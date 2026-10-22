#!/usr/bin/env python3
"""Find structural variants supported by both DELLY and SvABA.

Matching follows the same rules used by merge_delly_gridss.py:
  * both breakpoint distances are <= 300 bp (configurable);
  * both breakpoint strand orientations match;
  * intrachromosomal calls have >= 50% reciprocal overlap (configurable).

The shared-event TSV mirrors the DELLY/GRIDSS output layout, with SvABA
columns replacing GRIDSS columns.  SvABA ``tovcf --always-bnd`` output is
supported: intrachromosomal records retain the caller label BND, while
interchromosomal BND records are labelled TRA.  A second output VCF contains
the selected DELLY records, preserving the DELLY header, INFO, FORMAT, and
sample columns.

Only the Python standard library is required. BCF input additionally requires
bcftools for conversion to VCF text.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import re
import shutil
import subprocess
import sys
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


BND_RE = re.compile(r"([\[\]])([^:\[\]]+):(\d+)\1")
CT_RE = re.compile(r"^([35])to([35])$")


class InputFormatError(RuntimeError):
    """Raised when an input record cannot be represented by two breakpoints."""


@dataclass(frozen=True)
class Event:
    caller: str
    record_id: str
    chrom1: str
    pos1: int
    strand1: str
    chrom2: str
    pos2: int
    strand2: str
    svtype: str
    qual: str

    @property
    def intrachromosomal(self) -> bool:
        return self.chrom1 == self.chrom2


@dataclass(frozen=True)
class Candidate:
    delly_index: int
    gridss_index: int  # Kept for compatibility; this is the SvABA index here.
    distance1: int
    distance2: int
    reciprocal_overlap: Optional[float]

    @property
    def cost(self) -> int:
        overlap_penalty = 0
        if self.reciprocal_overlap is not None:
            overlap_penalty = int(round((1.0 - self.reciprocal_overlap) * 1000))
        return (
            max(self.distance1, self.distance2) * 1_000_000
            + (self.distance1 + self.distance2) * 1_000
            + overlap_penalty
        )


def warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def normalize_chrom(chrom: str) -> str:
    value = chrom.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    if value.upper() in {"M", "MT"}:
        value = "MT"
    return value.upper()


def chrom_sort_key(chrom: str) -> Tuple[int, object]:
    core = normalize_chrom(chrom)
    if core.isdigit():
        return 0, int(core)
    rank = {"X": 23, "Y": 24, "MT": 25}
    if core in rank:
        return 0, rank[core]
    return 1, core


def display_chrom(chrom: str) -> str:
    return "chr" + normalize_chrom(chrom)


def endpoint_key(chrom: str, pos: int) -> Tuple[Tuple[int, object], int]:
    return chrom_sort_key(chrom), pos


def normalize_svtype(svtype: str, interchromosomal: bool = False) -> str:
    value = (svtype or "NA").upper()
    if interchromosomal or value in {"TRA", "CTX", "ITX"}:
        return "TRA"
    # BND describes the VCF representation, not necessarily a translocation.
    # In particular, SvABA tovcf --always-bnd emits intrachromosomal deletions,
    # duplications and inversions as paired BND records as well.
    if value == "BND":
        return "BND"
    if value.startswith("DUP"):
        return "DUP"
    if value.startswith("DEL"):
        return "DEL"
    if value.startswith("INV"):
        return "INV"
    if value.startswith("INS"):
        return "INS"
    return value


def make_event(
    caller: str,
    record_id: str,
    chrom1: str,
    pos1: int,
    strand1: str,
    chrom2: str,
    pos2: int,
    strand2: str,
    svtype: str,
    qual: str,
) -> Event:
    left = (display_chrom(chrom1), int(pos1), strand1)
    right = (display_chrom(chrom2), int(pos2), strand2)
    if endpoint_key(right[0], right[1]) < endpoint_key(left[0], left[1]):
        left, right = right, left
    return Event(
        caller=caller,
        record_id=record_id,
        chrom1=left[0],
        pos1=left[1],
        strand1=left[2],
        chrom2=right[0],
        pos2=right[1],
        strand2=right[2],
        svtype=normalize_svtype(svtype, left[0] != right[0]),
        qual=qual,
    )


def parse_info(text: str) -> Dict[str, str]:
    info: Dict[str, str] = {}
    for item in text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = ""
    return info


def parse_ct(ct: str, svtype: str) -> Tuple[str, str]:
    match = CT_RE.match(ct or "")
    if match:
        return (
            "+" if match.group(1) == "3" else "-",
            "+" if match.group(2) == "3" else "-",
        )
    fallback = {"DEL": ("+", "-"), "DUP": ("-", "+"), "INS": ("+", "-")}
    value = normalize_svtype(svtype)
    if value in fallback:
        return fallback[value]
    raise InputFormatError(
        f"cannot determine orientation from CT={ct!r}, SVTYPE={svtype!r}"
    )


def parse_bnd_alt(alt: str) -> Tuple[str, int, str, str]:
    if "," in alt:
        raise InputFormatError(f"multi-allelic BND ALT is unsupported: {alt}")
    match = BND_RE.search(alt)
    if not match:
        raise InputFormatError(f"ALT is not a paired BND allele: {alt}")
    bracket, chrom2, pos2_text = match.groups()
    first_bracket = min(i for i in (alt.find("["), alt.find("]")) if i >= 0)
    local_strand = "+" if first_bracket > 0 else "-"
    remote_strand = "+" if bracket == "]" else "-"
    return chrom2, int(pos2_text), local_strand, remote_strand


@contextmanager
def open_variant_text(path: Path, bcftools: str) -> Iterator[Iterable[str]]:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if path.suffix.lower() == ".bcf":
        executable = shutil.which(bcftools)
        if executable is None:
            raise RuntimeError(
                f"bcftools executable {bcftools!r} was not found; activate the "
                "bcftools environment before reading BCF files"
            )
        process = subprocess.Popen(
            [executable, "view", "-Ov", "--no-version", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )
        assert process.stdout is not None
        assert process.stderr is not None
        try:
            yield process.stdout
        finally:
            process.stdout.close()
            stderr = process.stderr.read()
            return_code = process.wait()
            process.stderr.close()
            if return_code != 0:
                raise RuntimeError(
                    f"bcftools failed for {path} (exit {return_code}): {stderr.strip()}"
                )
    elif suffixes and suffixes[-1] == ".gz":
        handle = gzip.open(path, "rt", encoding="utf-8")
        try:
            yield handle
        finally:
            handle.close()
    else:
        handle = path.open("rt", encoding="utf-8")
        try:
            yield handle
        finally:
            handle.close()


def vcf_records(lines: Iterable[str]) -> Iterator[List[str]]:
    for line_number, line in enumerate(lines, 1):
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 8:
            raise InputFormatError(
                f"VCF line {line_number} has fewer than 8 columns"
            )
        yield fields


def is_pass(filter_value: str) -> bool:
    return filter_value == "PASS"


def mate_group_id(record_id: str, info: Dict[str, str]) -> str:
    mates = [value for value in info.get("MATEID", "").split(",") if value]
    if mates:
        return ",".join(sorted(set([record_id] + mates)))
    return record_id


def event_dedup_key(event: Event, group_id: str) -> Tuple[object, ...]:
    return (
        group_id,
        event.chrom1,
        event.pos1,
        event.strand1,
        event.chrom2,
        event.pos2,
        event.strand2,
    )


def reciprocal_overlap(first: Event, second: Event) -> float:
    first_start, first_end = sorted((first.pos1, first.pos2))
    second_start, second_end = sorted((second.pos1, second.pos2))
    if first_end == first_start:
        first_end += 1
    if second_end == second_start:
        second_end += 1
    overlap = max(0, min(first_end, second_end) - max(first_start, second_start))
    if overlap == 0:
        return 0.0
    return min(
        overlap / (first_end - first_start),
        overlap / (second_end - second_start),
    )


@dataclass
class FlowEdge:
    to: int
    reverse: int
    capacity: int
    cost: int
    candidate_index: Optional[int] = None


def add_flow_edge(
    graph: List[List[FlowEdge]],
    source: int,
    target: int,
    capacity: int,
    cost: int,
    candidate_index: Optional[int] = None,
) -> None:
    forward = FlowEdge(target, len(graph[target]), capacity, cost, candidate_index)
    reverse = FlowEdge(source, len(graph[source]), 0, -cost, None)
    graph[source].append(forward)
    graph[target].append(reverse)


def select_one_to_one(
    delly_count: int,
    svaba_count: int,
    candidates: Sequence[Candidate],
) -> set:
    """Return maximum-cardinality, minimum-cost bipartite matching."""
    if not candidates:
        return set()
    source = 0
    delly_offset = 1
    svaba_offset = delly_offset + delly_count
    sink = svaba_offset + svaba_count
    graph: List[List[FlowEdge]] = [[] for _ in range(sink + 1)]
    for index in {candidate.delly_index for candidate in candidates}:
        add_flow_edge(graph, source, delly_offset + index, 1, 0)
    for index in {candidate.gridss_index for candidate in candidates}:
        add_flow_edge(graph, svaba_offset + index, sink, 1, 0)
    for candidate_index, candidate in enumerate(candidates):
        add_flow_edge(
            graph,
            delly_offset + candidate.delly_index,
            svaba_offset + candidate.gridss_index,
            1,
            candidate.cost,
            candidate_index,
        )

    infinity = 10**30
    potential = [0] * len(graph)
    while True:
        distance = [infinity] * len(graph)
        previous_node = [-1] * len(graph)
        previous_edge = [-1] * len(graph)
        distance[source] = 0
        queue: List[Tuple[int, int]] = [(0, source)]
        while queue:
            current_distance, node = heapq.heappop(queue)
            if current_distance != distance[node]:
                continue
            for edge_index, edge in enumerate(graph[node]):
                if edge.capacity <= 0:
                    continue
                reduced = edge.cost + potential[node] - potential[edge.to]
                new_distance = current_distance + reduced
                if new_distance < distance[edge.to]:
                    distance[edge.to] = new_distance
                    previous_node[edge.to] = node
                    previous_edge[edge.to] = edge_index
                    heapq.heappush(queue, (new_distance, edge.to))
        if distance[sink] == infinity:
            break
        for node in range(len(graph)):
            if distance[node] < infinity:
                potential[node] += distance[node]
        node = sink
        while node != source:
            parent = previous_node[node]
            edge_index = previous_edge[node]
            edge = graph[parent][edge_index]
            edge.capacity -= 1
            graph[node][edge.reverse].capacity += 1
            node = parent

    selected = set()
    for index in {candidate.delly_index for candidate in candidates}:
        for edge in graph[delly_offset + index]:
            if edge.candidate_index is not None and edge.capacity == 0:
                selected.add(edge.candidate_index)
    return selected


def orientation_class(event: Event) -> str:
    if not event.intrachromosomal:
        return "TRA"
    mapping = {
        ("+", "-"): "DEL_LIKE",
        ("-", "+"): "DUP_LIKE",
        ("+", "+"): "INV_3TO3",
        ("-", "-"): "INV_5TO5",
    }
    return mapping[(event.strand1, event.strand2)]


OUTPUT_FIELDS = [
    "sample",
    "merge_id",
    "chrom1",
    "merged_pos1",
    "strand1",
    "chrom2",
    "merged_pos2",
    "strand2",
    "orientation_class",
    "intrachromosomal",
    "reciprocal_overlap",
    "breakpoint_distance1",
    "breakpoint_distance2",
    "delly_id",
    "delly_svtype",
    "delly_pos1",
    "delly_pos2",
    "delly_qual",
    "svaba_id",
    "svaba_svtype",
    "svaba_pos1",
    "svaba_pos2",
    "svaba_qual",
    "type_labels_agree",
]

CANDIDATE_FIELDS = OUTPUT_FIELDS + ["selected_one_to_one"]

SUMMARY_FIELDS = [
    "sample",
    "svaba_sv_count",
    "delly_sv_count",
    "shared_sv_count",
    "delly_records_seen",
    "delly_events_loaded",
    "svaba_records_seen",
    "svaba_events_loaded",
    "candidate_pairs",
    "shared_events",
    "unmatched_delly_events",
    "unmatched_svaba_events",
    "delly_skipped_filter",
    "svaba_skipped_filter",
    "delly_skipped_unsupported",
    "svaba_skipped_unsupported",
    "svaba_skipped_unknown_orientation",
    "delly_duplicate_mates",
    "svaba_duplicate_mates",
    "delly_vcf_records_written",
]


def svaba_group_id(record_id: str, info: Dict[str, str]) -> str:
    """Return one stable ID for a symbolic event or a pair of BND records."""
    return info.get("EVENT") or mate_group_id(record_id, info)


def load_delly_events_fixed(
    path: Path,
    bcftools: str,
    pass_only: bool,
    stats: Counter,
) -> List[Event]:
    """Load DELLY records, including CHR2+END translocations."""
    events: List[Event] = []
    seen: set = set()
    with open_variant_text(path, bcftools) as lines:
        for fields in vcf_records(lines):
            stats["delly_records_seen"] += 1
            chrom, pos_text, record_id, ref, alt, qual, filt, info_text = fields[:8]
            if pass_only and not is_pass(filt):
                stats["delly_skipped_filter"] += 1
                continue
            info = parse_info(info_text)
            svtype = info.get("SVTYPE", alt.strip("<>"))
            try:
                pos1 = int(pos_text)
                if normalize_svtype(svtype) in {"TRA", "BND"}:
                    chrom2 = info.get("CHR2")
                    pos2_text = info.get("POS2") or info.get("END")
                    if chrom2 and pos2_text not in {None, "", "."}:
                        pos2 = int(pos2_text)
                        alt_s1 = alt_s2 = ""
                    else:
                        chrom2, pos2, alt_s1, alt_s2 = parse_bnd_alt(alt)
                else:
                    chrom2 = chrom
                    end_text = info.get("END")
                    if end_text in {None, "", "."}:
                        pos2 = pos1 + max(1, len(ref) - 1)
                    else:
                        pos2 = int(end_text)
                    alt_s1 = alt_s2 = ""

                try:
                    strand1, strand2 = parse_ct(info.get("CT", ""), svtype)
                except InputFormatError:
                    if alt_s1 and alt_s2:
                        strand1, strand2 = alt_s1, alt_s2
                    else:
                        raise

                event = make_event(
                    "DELLY",
                    record_id,
                    chrom,
                    pos1,
                    strand1,
                    chrom2,
                    pos2,
                    strand2,
                    svtype,
                    qual,
                )
            except (ValueError, InputFormatError) as exc:
                stats["delly_skipped_unsupported"] += 1
                warn(f"DELLY {path.name} record {record_id}: {exc}; skipped")
                continue

            group_id = mate_group_id(record_id, info)
            key = event_dedup_key(event, group_id)
            if key in seen:
                stats["delly_duplicate_mates"] += 1
                continue
            seen.add(key)
            if group_id != record_id:
                event = Event(**{**event.__dict__, "record_id": group_id})
            events.append(event)

    stats["delly_events_loaded"] = len(events)
    return events


def load_svaba_events(
    path: Path,
    pass_only: bool,
    allow_symbolic_inv_compatible: bool,
    stats: Counter,
) -> List[Event]:
    """Load SvABA symbolic records or paired ``tovcf --always-bnd`` records."""
    events: List[Event] = []
    seen: set = set()
    with open_variant_text(path, bcftools="bcftools") as lines:
        for fields in vcf_records(lines):
            stats["svaba_records_seen"] += 1
            chrom, pos_text, record_id, ref, alt, qual, filt, info_text = fields[:8]
            if pass_only and not is_pass(filt):
                stats["svaba_skipped_filter"] += 1
                continue
            info = parse_info(info_text)
            svtype = normalize_svtype(info.get("SVTYPE", alt.strip("<>")))

            try:
                pos1 = int(pos_text)
                if alt.startswith("<") and alt.endswith(">"):
                    chrom2 = chrom
                    end_text = info.get("END")
                    if end_text in {None, "", "."}:
                        raise InputFormatError(f"symbolic {alt} record has no INFO/END")
                    pos2 = int(end_text)
                    if svtype == "DEL":
                        strand1, strand2 = "+", "-"
                    elif svtype == "DUP":
                        strand1, strand2 = "-", "+"
                    elif svtype == "INS":
                        strand1, strand2 = "+", "-"
                    elif svtype == "INV":
                        if not allow_symbolic_inv_compatible:
                            stats["svaba_skipped_unknown_orientation"] += 1
                            continue
                        # SvABA's symbolic <INV> omits whether the adjacency is
                        # 3-to-3 (++) or 5-to-5 (--). '?' is handled as a
                        # same-strand compatibility wildcard in find_candidates.
                        strand1 = strand2 = "?"
                    else:
                        raise InputFormatError(f"unsupported symbolic ALT {alt}")
                else:
                    chrom2, pos2, strand1, strand2 = parse_bnd_alt(alt)

                group_id = svaba_group_id(record_id, info)
                event = make_event(
                    "SVABA",
                    group_id,
                    chrom,
                    pos1,
                    strand1,
                    chrom2,
                    pos2,
                    strand2,
                    svtype,
                    qual,
                )
            except (ValueError, InputFormatError) as exc:
                stats["svaba_skipped_unsupported"] += 1
                warn(f"SvABA {path.name} record {record_id}: {exc}; skipped")
                continue

            key = event_dedup_key(event, group_id)
            if key in seen:
                stats["svaba_duplicate_mates"] += 1
                continue
            seen.add(key)
            events.append(event)

    stats["svaba_events_loaded"] = len(events)
    return events


def orientation_matches(delly: Event, svaba: Event) -> bool:
    if svaba.strand1 == "?" and svaba.strand2 == "?":
        return svaba.svtype == "INV" and delly.strand1 == delly.strand2
    return (
        delly.strand1 == svaba.strand1
        and delly.strand2 == svaba.strand2
    )


def find_candidates(
    delly: Sequence[Event],
    svaba: Sequence[Event],
    max_distance: int,
    min_overlap: float,
) -> List[Candidate]:
    """Find all pairs satisfying chromosome, orientation, distance and RO."""
    candidates: List[Candidate] = []
    for delly_index, d_event in enumerate(delly):
        for svaba_index, s_event in enumerate(svaba):
            if d_event.chrom1 != s_event.chrom1 or d_event.chrom2 != s_event.chrom2:
                continue
            if not orientation_matches(d_event, s_event):
                continue
            distance1 = abs(d_event.pos1 - s_event.pos1)
            distance2 = abs(d_event.pos2 - s_event.pos2)
            if distance1 > max_distance or distance2 > max_distance:
                continue
            ro: Optional[float] = None
            if d_event.intrachromosomal:
                ro = reciprocal_overlap(d_event, s_event)
                if ro + 1e-12 < min_overlap:
                    continue
            candidates.append(
                Candidate(delly_index, svaba_index, distance1, distance2, ro)
            )
    return candidates


def type_labels_agree(delly: Event, svaba: Event) -> bool:
    return delly.svtype == svaba.svtype


def candidate_row(
    sample: str,
    merge_number: int,
    candidate: Candidate,
    delly: Event,
    svaba: Event,
) -> Dict[str, object]:
    return {
        "sample": sample,
        "merge_id": f"{sample}_SHARED_{merge_number:06d}",
        "chrom1": delly.chrom1,
        "merged_pos1": (delly.pos1 + svaba.pos1) // 2,
        "strand1": delly.strand1,
        "chrom2": delly.chrom2,
        "merged_pos2": (delly.pos2 + svaba.pos2) // 2,
        "strand2": delly.strand2,
        "orientation_class": orientation_class(delly),
        "intrachromosomal": str(delly.intrachromosomal).lower(),
        "reciprocal_overlap": (
            "NA"
            if candidate.reciprocal_overlap is None
            else f"{candidate.reciprocal_overlap:.6f}"
        ),
        "breakpoint_distance1": candidate.distance1,
        "breakpoint_distance2": candidate.distance2,
        "delly_id": delly.record_id,
        "delly_svtype": delly.svtype,
        "delly_pos1": delly.pos1,
        "delly_pos2": delly.pos2,
        "delly_qual": delly.qual,
        "svaba_id": svaba.record_id,
        "svaba_svtype": svaba.svtype,
        "svaba_pos1": svaba.pos1,
        "svaba_pos2": svaba.pos2,
        "svaba_qual": svaba.qual,
        "type_labels_agree": str(type_labels_agree(delly, svaba)).lower(),
    }


def write_tsv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_selected_delly_vcf(
    source: Path,
    output: Path,
    selected_group_ids: set,
    bcftools: str,
    max_distance: int,
    min_overlap: float,
) -> int:
    """Write complete selected DELLY records using the original VCF schema."""
    output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    inserted_provenance = False
    with open_variant_text(source, bcftools) as lines, output.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        for line in lines:
            if line.startswith("#CHROM") and not inserted_provenance:
                handle.write(
                    "##merge_delly_svaba=<"
                    f"MaxBreakpointDistance={max_distance},"
                    f"MinReciprocalOverlap={min_overlap:g}>\n"
                )
                inserted_provenance = True
                handle.write(line)
                continue
            if line.startswith("#"):
                handle.write(line)
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            if mate_group_id(fields[2], info) in selected_group_ids:
                handle.write(line)
                written += 1
    return written


def merge_sample(
    sample: str,
    delly_path: Path,
    svaba_path: Path,
    output_tsv: Path,
    output_vcf: Path,
    bcftools: str,
    max_distance: int,
    min_overlap: float,
    pass_only: bool,
    allow_symbolic_inv_compatible: bool,
    candidates_path: Optional[Path] = None,
) -> Dict[str, object]:
    stats: Counter = Counter()
    delly_events = load_delly_events_fixed(delly_path, bcftools, pass_only, stats)
    svaba_events = load_svaba_events(
        svaba_path, pass_only, allow_symbolic_inv_compatible, stats
    )
    candidates = find_candidates(
        delly_events, svaba_events, max_distance, min_overlap
    )
    selected_indices = select_one_to_one(
        len(delly_events), len(svaba_events), candidates
    )
    selected_candidates = [candidates[index] for index in sorted(selected_indices)]
    selected_candidates.sort(
        key=lambda item: (
            chrom_sort_key(delly_events[item.delly_index].chrom1),
            delly_events[item.delly_index].pos1,
            chrom_sort_key(delly_events[item.delly_index].chrom2),
            delly_events[item.delly_index].pos2,
        )
    )

    selected_rows = [
        candidate_row(
            sample,
            number,
            candidate,
            delly_events[candidate.delly_index],
            svaba_events[candidate.gridss_index],
        )
        for number, candidate in enumerate(selected_candidates, 1)
    ]
    write_tsv(output_tsv, OUTPUT_FIELDS, selected_rows)

    if candidates_path is not None:
        candidate_rows = []
        for index, candidate in enumerate(candidates):
            row = candidate_row(
                sample,
                index + 1,
                candidate,
                delly_events[candidate.delly_index],
                svaba_events[candidate.gridss_index],
            )
            row["selected_one_to_one"] = str(index in selected_indices).lower()
            candidate_rows.append(row)
        write_tsv(candidates_path, CANDIDATE_FIELDS, candidate_rows)

    selected_delly_ids = {
        delly_events[candidate.delly_index].record_id
        for candidate in selected_candidates
    }
    stats["delly_vcf_records_written"] = write_selected_delly_vcf(
        delly_path,
        output_vcf,
        selected_delly_ids,
        bcftools,
        max_distance,
        min_overlap,
    )

    matched_delly = {candidate.delly_index for candidate in selected_candidates}
    matched_svaba = {candidate.gridss_index for candidate in selected_candidates}
    stats["candidate_pairs"] = len(candidates)
    stats["shared_events"] = len(selected_candidates)
    stats["unmatched_delly_events"] = len(delly_events) - len(matched_delly)
    stats["unmatched_svaba_events"] = len(svaba_events) - len(matched_svaba)

    summary = {field: stats.get(field, 0) for field in SUMMARY_FIELDS}
    summary["sample"] = sample
    summary["svaba_sv_count"] = len(svaba_events)
    summary["delly_sv_count"] = len(delly_events)
    summary["shared_sv_count"] = len(selected_candidates)
    return summary


def sample_from_name(path: Path, suffix: str, caller: str) -> str:
    if not path.name.endswith(suffix):
        raise ValueError(f"{path.name!r} does not end with {caller} suffix {suffix!r}")
    return path.name[: -len(suffix)]


def add_matching_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--bcftools", default="bcftools", help="bcftools executable (default: bcftools)"
    )
    parser.add_argument(
        "--max-breakpoint-distance",
        type=int,
        default=300,
        metavar="BP",
        help="maximum distance at each breakpoint (default: 300)",
    )
    parser.add_argument(
        "--min-reciprocal-overlap",
        type=float,
        default=0.5,
        metavar="FRACTION",
        help="minimum overlap fraction for each intrachromosomal interval (default: 0.5)",
    )
    parser.add_argument(
        "--allow-nonpass",
        action="store_true",
        help="include non-PASS records (default: only FILTER=PASS records)",
    )
    parser.add_argument(
        "--allow-symbolic-inv-compatible",
        action="store_true",
        help=(
            "allow SvABA symbolic <INV> to match either DELLY same-strand orientation; "
            "default excludes it because the exact ++/-- orientation is absent"
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Select DELLY/SvABA structural variants supported by both callers."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    single = subparsers.add_parser("single", help="merge one DELLY/SvABA sample pair")
    single.add_argument("--sample", help="sample label (default: inferred from DELLY filename)")
    single.add_argument("--delly", type=Path, required=True, help="DELLY somatic VCF or BCF")
    single.add_argument("--svaba", type=Path, required=True, help="filtered SvABA SV VCF")
    single.add_argument("--output", type=Path, required=True, help="shared-event TSV")
    single.add_argument(
        "--output-vcf", type=Path, help="selected DELLY-format VCF (default: beside TSV)"
    )
    single.add_argument("--candidates-output", type=Path, help="optional all-candidate QC TSV")
    single.add_argument("--delly-suffix", default=".somatic.vcf")
    add_matching_options(single)

    batch = subparsers.add_parser("batch", help="merge all paired files in two directories")
    batch.add_argument("--delly-dir", type=Path, required=True)
    batch.add_argument("--svaba-dir", type=Path, required=True)
    batch.add_argument("--output-dir", type=Path, required=True)
    batch.add_argument("--delly-suffix", default=".somatic.vcf")
    batch.add_argument("--svaba-suffix", default=".somatic.vaf0.20.ad4.sv.vcf")
    batch.add_argument(
        "--write-candidates", action="store_true", help="write an all-candidate QC TSV per sample"
    )
    add_matching_options(batch)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.max_breakpoint_distance < 0:
        raise ValueError("--max-breakpoint-distance must be >= 0")
    if not 0 <= args.min_reciprocal_overlap <= 1:
        raise ValueError("--min-reciprocal-overlap must be between 0 and 1")


def run_single(args: argparse.Namespace) -> int:
    if not args.delly.is_file():
        raise FileNotFoundError(args.delly)
    if not args.svaba.is_file():
        raise FileNotFoundError(args.svaba)
    sample = args.sample or (
        sample_from_name(args.delly, args.delly_suffix, "DELLY")
        if args.delly.name.endswith(args.delly_suffix)
        else args.delly.stem
    )
    output_vcf = args.output_vcf or args.output.with_name(
        args.output.name.removesuffix(".tsv") + ".delly.vcf"
    )
    summary = merge_sample(
        sample,
        args.delly,
        args.svaba,
        args.output,
        output_vcf,
        args.bcftools,
        args.max_breakpoint_distance,
        args.min_reciprocal_overlap,
        not args.allow_nonpass,
        args.allow_symbolic_inv_compatible,
        args.candidates_output,
    )
    summary_path = args.output.with_name(args.output.stem + ".summary.tsv")
    write_tsv(summary_path, SUMMARY_FIELDS, [summary])
    print(
        f"{sample}: {summary['shared_sv_count']} shared events; "
        f"wrote {args.output}, {output_vcf}, and {summary_path}",
        file=sys.stderr,
    )
    return 0


def run_batch(args: argparse.Namespace) -> int:
    if not args.delly_dir.is_dir():
        raise NotADirectoryError(args.delly_dir)
    if not args.svaba_dir.is_dir():
        raise NotADirectoryError(args.svaba_dir)
    delly_files = sorted(args.delly_dir.glob("*" + args.delly_suffix))
    if not delly_files:
        raise FileNotFoundError(
            f"no files ending in {args.delly_suffix!r} under {args.delly_dir}"
        )

    summaries = []
    for delly_path in delly_files:
        sample = sample_from_name(delly_path, args.delly_suffix, "DELLY")
        svaba_path = args.svaba_dir / (sample + args.svaba_suffix)
        if not svaba_path.is_file():
            warn(f"{sample}: SvABA file not found ({svaba_path}); skipped")
            continue

        output_tsv = args.output_dir / f"{sample}.delly_svaba.shared.tsv"
        output_vcf = args.output_dir / f"{sample}.delly_svaba.shared.delly.vcf"
        candidates_path = (
            args.output_dir / f"{sample}.delly_svaba.candidates.tsv"
            if args.write_candidates
            else None
        )
        summary = merge_sample(
            sample,
            delly_path,
            svaba_path,
            output_tsv,
            output_vcf,
            args.bcftools,
            args.max_breakpoint_distance,
            args.min_reciprocal_overlap,
            not args.allow_nonpass,
            args.allow_symbolic_inv_compatible,
            candidates_path,
        )
        summaries.append(summary)
        print(f"{sample}: {summary['shared_sv_count']} shared events", file=sys.stderr)

    if not summaries:
        raise RuntimeError("no DELLY/SvABA sample pairs were found")

    total: Dict[str, object] = {field: 0 for field in SUMMARY_FIELDS}
    total["sample"] = "TOTAL"
    for summary in summaries:
        for field in SUMMARY_FIELDS:
            if field != "sample":
                total[field] = int(total[field]) + int(summary[field])
    summaries.append(total)
    summary_path = args.output_dir / "delly_svaba.merge_summary.tsv"
    write_tsv(summary_path, SUMMARY_FIELDS, summaries)
    print(
        f"TOTAL: SvABA={total['svaba_sv_count']}, "
        f"DELLY={total['delly_sv_count']}, shared={total['shared_sv_count']}",
        file=sys.stderr,
    )
    print(f"Wrote batch summary: {summary_path}", file=sys.stderr)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        validate_args(args)
        return run_single(args) if args.command == "single" else run_batch(args)
    except (
        FileNotFoundError,
        NotADirectoryError,
        InputFormatError,
        RuntimeError,
        ValueError,
    ) as exc:
        parser.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
