use strict;
use warnings;
use utf8;
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw( mkpath );
use lib $Bin;

my ($snp , $indel , $help);
my $T_depth_snv =5;
my $T_depth_Var_snv =2;
my $N_depth_snv = 5;
my $N_depth_Var_snv = 2;
my $T_Var_Strand_snv = 2;
my $T_Var_Frac_snv =0.04;
my $N_Var_Frac_snv =0.03;
my $somaticp_snv =0.05;
my $qss_vaf_cutoff_snv = 1.3;


my $T_depth_indel =5;
my $T_depth_Var_indel =2;
my $N_depth_indel = 5;
my $N_depth_Var_indel = 2;
my $T_Var_Strand_indel = 2;
my $T_Var_Frac_indel =0.04;
my $N_Var_Frac_indel =0.03;
my $somaticp_indel =0.05;
my $qss_vaf_cutoff_indel = 1.3;

GetOptions(
	"s|snp=s"	=>	\$snp,
	"i|indel=s"	=>	\$indel,
        "total_depth_snv=i"   =>    \$T_depth_snv,
        "var_depth_snv=i"     =>      \$T_depth_Var_snv,
        "vaf_snv=f"     =>      \$T_Var_Frac_snv,
        "total_depth_indel=i" =>     \$T_depth_indel,
        "var_depth_indel=i"    =>      \$T_depth_Var_indel,
        "vaf_indel=f"     =>     \$T_Var_Frac_indel,
	"help"	=>	\$help,
);

#if ($help or ! $snp){
#	&help;
#	exit;
#}

my $fn = 0;
my $file;

if(defined($snp) and !(defined($indel))){
	$fn = 1;
	$file = $snp;
}
elsif(!(defined($snp)) and defined($indel)){
	$fn=2;
	$file = $indel;
}
elsif(defined($snp) and defined($indel)){
	$fn=3;
}

my %pv = ('A'=>4,'C'=>5,'G'=>6,'T'=>7);  
if($fn ==1 or $fn == 2){
	if ($file =~ /\.gz/){
		open IN , "zcat $file|"; 
	}else{
		open IN , "$file"; 
	}
	while (<IN>){
		if (/^##/){
			print $_; 
			next;
		}elsif (/^#/){
			print '#' , join("\t" , qw/CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT	NORMAL  TUMOR/) , "\n";
			next;
		}
		chomp;
		my @F = split /\t/ , $_; 
		next unless $F[6] eq 'PASS';
		#next if exists $pon{"$F[0]\t$F[1]"};
		my @param=split(/;/,$F[7]); 
		my $p=0;
		my $nt = '';
		my $sgt = '';
		my $qss=0;
		while($p<=$#param){ 
			my @info=split(/[=,]/,$param[$p]);
			if ($info[0] eq 'NT'){ 
				$nt = $info[1];
			}elsif ($info[0] eq 'SGT'){
				$sgt = $info[1];
			}elsif ($info[0] eq 'QSS'){ 
				$qss = $info[1];
			}elsif($info[0] eq 'QSI'){
				$qss = $info[1];
			}
			$p++;
		}
		next unless $nt eq 'ref';
		my ($sr , $sa) = split /\-\>/ , $sgt;
		next if $sr eq $sa;
		my ($sa1 , $sa2) = split // , $sa;
		my $gt = '';
		if ($sa1 eq $sa2){
			$gt = '1/1';
		}else{
			$gt = '0/1';
		}
		#DP:FDP:SDP:SUBDP:AU:CU:GU:TU    60:2:0:0:0,0:0,0:0,0:58,60      294:13:0:0:0,0:2,2:11,16:268,279
		#DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50	22:22:22,23:0,1:0,16:19.25:2.20:0.00:0.11	51:51:36,41:9,15:12,88:51.07:9.64:0.00:0.18

		my @n=split(/:/,$F[9]);
		my @t=split(/:/,$F[10]);
		my ($nvar,$nref,$tvar,$tref,$ndp,$tdp);
		if ($fn==1){
			$nvar = $n[$pv{$F[4]}]; 
			$nref = $n[$pv{$F[3]}];
			$tvar = $t[$pv{$F[4]}];
			$tref = $t[$pv{$F[3]}];
		}else{
			$nvar = $n[3];
			$nref = $n[2];
			$tvar = $t[3];
			$tref = $t[2];
		}
		$nvar =~ s/,\d+//;
		$nref =~ s/,\d+//;
		$tvar =~ s/,\d+//;
		$tref =~ s/,\d+//;
		$ndp = $nvar + $nref;
		$tdp = $tvar + $tref;
		next if $ndp <= 1 or $tdp <= 1;

		my $tvf = sprintf("%.4f" , $tvar/$tdp);
		my $nvf = sprintf("%.4f" , $nvar/$ndp);
		my $qss_vaf;
		$qss_vaf = $qss * $tvf;
		#print($qss_vaf);

		$F[8] = "GT:AD:AF:$F[8]:QSSVAF"; 
		$F[9] = "0/0:$nref,$nvar:$nvf:$F[9]";
		$F[10] = "$gt:$tref,$tvar:$tvf:$F[10]:$qss_vaf"; 
		

		if($fn == 1 and $tdp>=$T_depth_snv and $tvar>=$T_depth_Var_snv and $ndp>=$N_depth_snv and $nvf<$N_Var_Frac_snv and $tvf>=$T_Var_Frac_snv and $nvf<$tvf/5 and $qss_vaf > $qss_vaf_cutoff_snv){ # 筛选出tumor深度大于8, 肿瘤突变覆盖深度大于4，正常样本深度大于8，正常的突变比例小于0.01，肿瘤的突变比例大于0.03，正常的vaf小于肿瘤的vaf的1/5
			print join("\t" , @F) , "\n";
		}
		if($fn == 2 and $tdp>=$T_depth_indel and $tvar>=$T_depth_Var_indel and $ndp>=$N_depth_indel and $nvf<$N_Var_Frac_indel and $tvf>=$T_Var_Frac_indel and $nvf<$tvf/5 and $qss_vaf > $qss_vaf_cutoff_indel){
			print join("\t" , @F) , "\n";
		}
	}
	close IN;
}
elsif($fn == 3){
	for my $file ($snp , $indel){ 
	$fn++;
	if ($file =~ /\.gz/){
		open IN , "zcat $file|";
	}else{
		open IN , "$file"; 
	}
	while (<IN>){
		if (/^##/){
			print $_ if $fn == 4; 
			next;
		}elsif (/^#/){
			print '#' , join("\t" , qw/CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT	NORMAL  TUMOR/) , "\n" if $fn == 4;
			next;
		}
		chomp;
		my @F = split /\t/ , $_; 
		next unless $F[6] eq 'PASS';
		#next if exists $pon{"$F[0]\t$F[1]"};
		my @param=split(/;/,$F[7]); 
		my $p=0;
		my $nt = '';
		my $sgt = '';
		my $qss=0;
		while($p<=$#param){
			my @info=split(/[=,]/,$param[$p]);
			if ($info[0] eq 'NT'){ 
				$nt = $info[1];
			}elsif ($info[0] eq 'SGT'){
				$sgt = $info[1];
			}elsif ($info[0] eq 'QSS'){
				$qss = $info[1];
			}elsif($info[0] eq 'QSI'){
				$qss = $info[1];
			}
			$p++;
		}
		next unless $nt eq 'ref';
		my ($sr , $sa) = split /\-\>/ , $sgt;
		next if $sr eq $sa;
		my ($sa1 , $sa2) = split // , $sa;
		my $gt = '';
		if ($sa1 eq $sa2){
			$gt = '1/1';
		}else{
			$gt = '0/1';
		}

		my @n=split(/:/,$F[9]);
		my @t=split(/:/,$F[10]);
		my ($nvar,$nref,$tvar,$tref,$ndp,$tdp);
		if ($fn==4){ 
			$nvar = $n[$pv{$F[4]}]; 
			$nref = $n[$pv{$F[3]}]; 
			$tvar = $t[$pv{$F[4]}]; 
			$tref = $t[$pv{$F[3]}];
		}else{
			$nvar = $n[3];
			$nref = $n[2]; 
			$tvar = $t[3];
			$tref = $t[2];
		}
		$nvar =~ s/,\d+//;
		$nref =~ s/,\d+//;
		$tvar =~ s/,\d+//;
		$tref =~ s/,\d+//;
		$ndp = $nvar + $nref;
		$tdp = $tvar + $tref;
		next if $ndp <= 1 or $tdp <= 1;

		my $tvf = sprintf("%.4f" , $tvar/$tdp);
		my $nvf = sprintf("%.4f" , $nvar/$ndp);
		my $qss_vaf;
		$qss_vaf = $qss * $tvf;

		$F[8] = "GT:AD:AF:$F[8]:QSSVAF";
		$F[9] = "0/0:$nref,$nvar:$nvf:$F[9]";
		$F[10] = "$gt:$tref,$tvar:$tvf:$F[10]:$qss_vaf"; # AD是
		

		if($fn == 4 and $tdp>=$T_depth_snv and $tvar>=$T_depth_Var_snv and $ndp>=$N_depth_snv and $nvf<$N_Var_Frac_snv and $tvf>=$T_Var_Frac_snv and $nvf<$tvf/5 and $qss_vaf > $qss_vaf_cutoff_snv){ # 筛选出tumor深度大于8, 肿瘤突变覆盖深度大于4，正常样本深度大于8，正常的突变比例小于0.01，肿瘤的突变比例大于0.03，正常的vaf小于肿瘤的vaf的1/5
			print join("\t" , @F) , "\n";
		}
		elsif($fn == 5 and $tdp>=$T_depth_indel and $tvar>=$T_depth_Var_indel and $ndp>=$N_depth_indel and $nvf<$N_Var_Frac_indel and $tvf>=$T_Var_Frac_indel and $nvf<$tvf/5 and $qss_vaf > $qss_vaf_cutoff_indel){
			print join("\t" , @F) , "\n";
		}
	}
	close IN;
}
}


sub help{
print << "EOF!";

EOF!
}
