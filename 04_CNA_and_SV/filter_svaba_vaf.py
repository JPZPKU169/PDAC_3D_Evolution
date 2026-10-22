#!/usr/bin/env python3
"""Filter a SvABA SV VCF by tumor VAF, tumor AD, and INFO/SOMATIC.

SvABA reports one AD and one DP value per sample.  Its allele fraction is
AD / DP (capped at 1 for SV records, where AD can occasionally exceed DP).
Paired BND records are selected as a complete EVENT so that one mate is not
left behind.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


@dataclass
class VcfRecord:
    line: str
    event_key: str
    has_somatic: bool
    tumor_ad: int | None
    tumor_dp: int | None
    tumor_vaf: float | None
    passes: bool


def open_text(path: Path, mode: str) -> TextIO:
    """Open plain text or gzip-compressed input/output by filename suffix."""
    if path.name.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter a SvABA VCF for INFO/SOMATIC, tumor AD, and tumor VAF "
            "(VAF=min(AD/DP,1)). Paired BND records sharing EVENT are retained together."
        )
    )
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input .vcf or .vcf.gz")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output .vcf or .vcf.gz")
    parser.add_argument(
        "-t",
        "--tumor-sample",
        help="Tumor sample name. If omitted, infer the unique sample whose name starts with 't'.",
    )
    parser.add_argument("--min-vaf", type=float, default=0.20, help="Minimum tumor VAF (default: 0.20)")
    parser.add_argument("--min-ad", type=int, default=4, help="Minimum tumor AD (default: 4)")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing output file")
    args = parser.parse_args()

    if not 0 <= args.min_vaf <= 1:
        parser.error("--min-vaf must be between 0 and 1")
    if args.min_ad < 0:
        parser.error("--min-ad must be non-negative")
    return args


def parse_info(text: str) -> dict[str, str | None]:
    info: dict[str, str | None] = {}
    if text == ".":
        return info
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = None
    return info


def parse_nonnegative_int(value: str | None) -> int | None:
    if value in (None, "", "."):
        return None
    try:
        number = int(value)
    except ValueError:
        return None
    return number if number >= 0 else None


def calculate_vaf(ad: int | None, dp: int | None) -> float | None:
    if ad is None or dp is None:
        return None
    if dp == 0:
        return 1.0 if ad > 0 else None
    return min(ad / dp, 1.0)


def make_event_key(fields: list[str], info: dict[str, str | None]) -> str:
    event = info.get("EVENT")
    if event:
        return "EVENT:" + event

    record_id = fields[2]
    mate_ids = info.get("MATEID")
    if record_id != "." and mate_ids:
        members = sorted({record_id, *mate_ids.split(",")})
        return "MATES:" + ",".join(members)
    if record_id != ".":
        return "ID:" + record_id
    return "RECORD:" + ":".join((fields[0], fields[1], fields[3], fields[4]))


def infer_tumor_sample(samples: list[str], requested: str | None) -> str:
    if requested:
        if requested not in samples:
            raise ValueError(
                f"Tumor sample {requested!r} is not in the VCF. Available samples: {', '.join(samples)}"
            )
        return requested

    candidates = [sample for sample in samples if sample.lower().startswith("t")]
    if len(candidates) == 1:
        return candidates[0]
    if len(samples) == 1:
        return samples[0]
    raise ValueError(
        "Could not infer one tumor sample. Use --tumor-sample. "
        f"Available samples: {', '.join(samples)}"
    )


def load_vcf(
    input_path: Path,
    requested_tumor: str | None,
    min_vaf: float,
    min_ad: int,
) -> tuple[list[str], list[VcfRecord], str]:
    headers: list[str] = []
    records: list[VcfRecord] = []
    samples: list[str] | None = None
    tumor_sample: str | None = None
    tumor_column: int | None = None

    with open_text(input_path, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##"):
                headers.append(line)
                continue
            if line.startswith("#CHROM"):
                headers.append(line)
                columns = line.rstrip("\r\n").split("\t")
                samples = columns[9:]
                if not samples:
                    raise ValueError("The VCF has no sample columns")
                tumor_sample = infer_tumor_sample(samples, requested_tumor)
                tumor_column = 9 + samples.index(tumor_sample)
                continue
            if line.startswith("#"):
                headers.append(line)
                continue
            if not line.strip():
                continue
            if tumor_column is None or tumor_sample is None:
                raise ValueError("A data record was found before the #CHROM header")

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) <= tumor_column:
                raise ValueError(f"Malformed VCF record at line {line_number}")

            info = parse_info(fields[7])
            format_keys = fields[8].split(":")
            sample_values = fields[tumor_column].split(":")
            sample_data = dict(zip(format_keys, sample_values))
            ad = parse_nonnegative_int(sample_data.get("AD"))
            dp = parse_nonnegative_int(sample_data.get("DP"))
            vaf = calculate_vaf(ad, dp)
            has_somatic = "SOMATIC" in info
            passes = (
                has_somatic
                and ad is not None
                and ad >= min_ad
                and vaf is not None
                and vaf >= min_vaf
            )
            records.append(
                VcfRecord(
                    line=line,
                    event_key=make_event_key(fields, info),
                    has_somatic=has_somatic,
                    tumor_ad=ad,
                    tumor_dp=dp,
                    tumor_vaf=vaf,
                    passes=passes,
                )
            )

    if samples is None or tumor_sample is None:
        raise ValueError("The VCF is missing a #CHROM header")
    return headers, records, tumor_sample


def write_filtered_vcf(
    output_path: Path,
    headers: list[str],
    records: list[VcfRecord],
    tumor_sample: str,
    min_vaf: float,
    min_ad: int,
) -> tuple[int, int, int]:
    selected_events = {record.event_key for record in records if record.passes}
    event_pass_states: dict[str, set[bool]] = {}
    for record in records:
        event_pass_states.setdefault(record.event_key, set()).add(record.passes)
    partial_events = sum(states == {False, True} for states in event_pass_states.values())

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(output_path, "w") as handle:
        provenance = (
            "##svaba_vaf_filter=<"
            f"TumorSample={tumor_sample},MinTumorVAF={min_vaf:g},"
            f"MinTumorAD={min_ad},RequireSomatic=true>\n"
        )
        inserted = False
        for header in headers:
            if header.startswith("#CHROM") and not inserted:
                handle.write(provenance)
                inserted = True
            handle.write(header)
        for record in records:
            if record.event_key in selected_events:
                handle.write(record.line)

    selected_records = sum(record.event_key in selected_events for record in records)
    return len(selected_events), selected_records, partial_events


def main() -> int:
    args = parse_args()
    input_path = args.input.resolve()
    output_path = args.output.resolve()

    if not input_path.is_file():
        print(f"ERROR: input file does not exist: {input_path}", file=sys.stderr)
        return 2
    if input_path == output_path:
        print("ERROR: input and output paths must be different", file=sys.stderr)
        return 2
    if output_path.exists() and not args.force:
        print(f"ERROR: output already exists (use --force): {output_path}", file=sys.stderr)
        return 2

    try:
        headers, records, tumor_sample = load_vcf(
            input_path,
            args.tumor_sample,
            args.min_vaf,
            args.min_ad,
        )
        selected_events, selected_records, partial_events = write_filtered_vcf(
            output_path,
            headers,
            records,
            tumor_sample,
            args.min_vaf,
            args.min_ad,
        )
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    all_events = {record.event_key for record in records}
    somatic_events = {record.event_key for record in records if record.has_somatic}
    print(f"Input: {input_path}")
    print(f"Output: {output_path}")
    print(f"Tumor sample: {tumor_sample}")
    print(f"Criteria: SOMATIC; tumor VAF >= {args.min_vaf:g}; tumor AD >= {args.min_ad}")
    print(f"Input records: {len(records)}")
    print(f"Input SV events: {len(all_events)}")
    print(f"SOMATIC events: {len(somatic_events)}")
    print(f"Selected records: {selected_records}")
    print(f"Selected SV events: {selected_events}")
    if partial_events:
        print(
            f"Note: {partial_events} paired event(s) had discordant per-record filter states; "
            "all mate records were retained.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
