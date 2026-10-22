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
my $T_depth =5;
my $T_depth_Var =2;
my $N_depth =5;
my $N_depth_Var =2;
my $T_Var_Strand =2;
my $T_Var_Frac =0.04;
my $N_Var_Frac =0.03;
my $somaticp =0.05;



my $nor_index = 9;
my $tum_index = 10;
GetOptions(
	"s|snp=s"	=>	\$snp,
        "total_depth=i"     =>    \$T_depth,
        "var_depth=i"  =>    \$T_depth_Var,
	"vaf=f"   =>      \$T_Var_Frac,
	"help"	=>	\$help,
);

if ($help or ! $snp){
	&help;
	exit;
}

my %pv = ('A'=>4,'C'=>5,'G'=>6,'T'=>7);
my $file = $snp;
if ($file =~ /\.gz/){
		open IN , "zcat $file|"; 
	}else{
		open IN , "$file"; 
	}
while (<IN>){
		if (/^##/){
			print $_ ; 
			next;
		}elsif (/^#/){
			my @file_head=split /\t/ , $_;
			if ($file_head[10] =~ /NOR/){  
				$nor_index = 10;
				$tum_index = 9;
			}
			#print $nor_index, $tum_index;
			print '#' , join("\t" , qw/CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT    NORMAL  TUMOR/) , "\n" ;
			next;
		}
		chomp;
		my @F = split /\t/ , $_; 
		next if length($F[3])!=length($F[4]);
        
		#next unless $F[6] eq 'PASS';
		#next if exists $pon{"$F[0]\t$F[1]"};
		
		# FORMAT: GT:AD:AF:DP:F1R2:F2R1:FAD:SB
		# 
		my @n=split(/:/,$F[$nor_index]);
		my @t=split(/:/,$F[$tum_index]);
		
		
		my ($nvar,$nref,$tvar,$tref);
        my ($ndp, $tdp);
		my ($nvaf, $tvaf);
        #$nvaf = $n[2];
        #$tvaf = $t[2];
        
        my @tdp_list = split(/,/,$t[1]);
		$tvar = $tdp_list[1];
        $tref = $tdp_list[0];
        $tdp = $tref + $tvar;
		$tvaf = sprintf("%.4f" , $tvar / $tdp);
        my @ndp_list = split(/,/,$n[1]);
		$nvar = $ndp_list[1];
		$nref = $ndp_list[0];
        $ndp = $nref + $nvar;
		$nvaf = sprintf("%.4f" , $nvar / $ndp);

        

		$F[9] = "$n[0]:$nref,$nvar:$nvaf:$ndp:$n[4]:$n[5]:$n[6]:$n[7]";
		$F[10] = "$t[0]:$tref,$tvar:$tvaf:$tdp:$t[4]:$t[5]:$t[6]:$t[7]";
        
        if($tdp >=$T_depth and $ndp >= $N_depth and $nvaf <$N_Var_Frac and $nvar <= $N_depth_Var and $tvaf >= $T_Var_Frac and $tvar >= $T_depth_Var){
			if(length($F[3]) > 1 and length($F[4]) > 1){ 
		    		for(my $loc = 0; $loc < length($F[3]); $loc = $loc + 1){
					my @F_tmp = @F;
					$F_tmp[1] = $F[1] + $loc;
		    		$F_tmp[3] = substr($F[3], $loc, 1);
					$F_tmp[4] = substr($F[4], $loc, 1);
					print join("\t", @F_tmp) , "\n";
				}
			}
			else{
				print join("\t", @F) , "\n";
			}
			#print join("\t", @F) , "\n";
        }

	}
close IN;
