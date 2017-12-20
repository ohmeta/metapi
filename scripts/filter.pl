#!/usr/bin/perl -w
use warnings;
use strict;
use File::Basename; 

die &usage if @ARGV < 6;
sub usage {
	print <<USAGE;
usage:
se pattern:
	perl $0 fq1 <prefix> <Qual system(33|64)> <min length> <seed OA> <fragment OA>
e.g	perl $0 sample.fq clean 33 30 0.99 0.90
	perl $0 sample.1.fq,sample.2.fq clean 33 30 0.99 0.90
USAGE
	exit;
}
sub openMethod {shift;return(($_=~/\.gz$/)?"gzip -dc $_|":"$_")}

### BODY ###
my ($fq,$pfx,$Qsys,$minLen,$Scut,$Qcut) = @ARGV;
$Qsys ||= 33;
my @fqs = split /,/,$fq;
my $mode= (@fqs==2)?"PE":"SE";

open F1,"gzip -dc $fqs[0] |",or die "error\n";
if($mode eq "PE"){
	open F1,"gzip -dc $fqs[0] |",or die "error\n";
	open O1,"|gzip >$pfx.clean.1.fq.gz" or die "error OUT\n";
	open F2,"gzip -dc $fqs[1] |",or die "error\n";
	open O2,"|gzip >$pfx.clean.2.fq.gz" or die "error OUT\n";
	open O3,"|gzip >$pfx.clean.single.fq.gz" or die "error OUT\n";
}else{
	open F1,"gzip -dc $fqs[0] |",or die "error\n";
	open O1,"|gzip >$pfx.clean.fq.gz" or die "error OUT\n";
}
open STAT,"> $pfx.clean.stat_out",or die "error\n";

my %STAT;
my ($total, @remainQ, @sum_bp, @sum_oa, @sum_Q, @sum_s) = ();
my @min_bp = (1e9,1e9,1e9);
my @max_bp = (0, 0, 0);

my (@fqID,@seq,@num,@qual,@originLength,@Tlength,@Aqual,@PQ,@start,@len,@count,@avgQ) =();
while(<F1>){
	#F1 info
	(@fqID,@seq,@num,@qual,@originLength,@Tlength,@Aqual,@PQ,@start,@len,@count,@avgQ) =();
	chomp;
	($fqID[0],$seq[0],$num[0],$qual[0])= &fqRead($_,\*F1);
	if($mode eq "SE"){
		my @a = split /\t| /,$fqID[0];
		($fqID[0] = $a[0]) =~ s/\/[12]$//;
	}else{
		$total ++;
		my $l2 = <F2>;
		($fqID[1],$seq[1],$num[1],$qual[1])= &fqRead($l2,\*F2);
	}
	$total ++;
	
	# trim & filter
	if($mode eq "PE"){
		($Aqual[0],$PQ[0],$start[0],$len[0],$seq[0],$qual[0]) = &ca1_cut($Scut,$Qcut,$Qsys,$seq[0],$qual[0]);
		($Aqual[1],$PQ[1],$start[1],$len[1],$seq[1],$qual[1]) = &ca1_cut($Scut,$Qcut,$Qsys,$seq[1],$qual[1]);
		 # filter
		if( $len[0] >= $minLen && $Aqual[0] >= $Qcut) {
			if( $len[1] >= $minLen && $Aqual[1] >= $Qcut) {
				print O1 "$fqID[0] length=$len[0]\n$seq[0]\n$num[0]\n$qual[0]\n";
				print O2 "$fqID[1] length=$len[1]\n$seq[1]\n$num[1]\n$qual[1]\n";
				&cumulate(0,0);
				&cumulate(1,1);
			}else{
				print O3 "$fqID[0] length=$len[0]\n$seq[0]\n$num[0]\n$qual[0]\n";
				&cumulate(2,0);
			}
		}else{
			if( $len[1] >= $minLen && $Aqual[1] >= $Qcut) {
				print O3 "$fqID[1] length=$len[1]\n$seq[1]\n$num[1]\n$qual[1]\n";
				&cumulate(2,0);
			}
		}
	}else{
		($Aqual[0],$PQ[0],$start[0],$len[0],$seq[0],$qual[0]) = &ca1_cut($Scut,$Qcut,$Qsys,$seq[0],$qual[0]);
	    # filter
	    if( $len[0] >= $minLen && $Aqual[0] >= $Qcut) {
        	print O1 "$fqID[0] length=$len[0]\n$seq[0]\n$num[0]\n$qual[0]\n";
	        # stat
			&cumulate(0,0);
	    }
	}
}
close F1;
close O1;
if($mode eq "PE"){
	close F2;
	close O2;
	close O3;
}
my ($max_bp, $min_bp, $sum_bp, $sum_oa, $sum_Q, $sum_s, $remainQ) = (0, 1e9, 0, 0, 0, 0, 0);
for(my $i=0;$i<@remainQ;$i++){
	$max_bp  = ($max_bp > $max_bp[$i])?$max_bp:$max_bp[$i];
	$min_bp  = ($min_bp < $min_bp[$i])?$min_bp:$min_bp[$i];
	$sum_bp += $sum_bp[$i];
	$sum_oa += $sum_oa[$i];
	$sum_Q  += $sum_Q[$i];
	$sum_s  += $sum_s[$i];
	$remainQ+= $remainQ[$i];
}
my $avgL = sprintf("%.0f",$sum_bp / $remainQ);
my $avgOA= sprintf("%.4f",$sum_oa / $remainQ);
my $avgQ = sprintf("%.0f",$sum_Q  / $sum_bp);
my $avgS = sprintf("%.0f",$sum_s  / $remainQ);
my $rate = sprintf("%.4f",$remainQ / $total);
my $tag = basename($pfx);
#my $debugHead = ($debug)?"\tN>$n|Len<$lf|PQ<$Qf|N+Len|N+PQ|Len+PQ|HOMER":"";

print STAT "Total\tmax\tmin\tavgLen\tavgStart\tavgOA\tavgQ\tremain\trate\n";
print STAT "$total\t$max_bp\t$min_bp\t$avgL\t$avgS\t$avgOA\t$avgQ\t$remainQ\t$rate\n";
if($mode eq "PE"){
	
	print STAT "Read_1\tmax\tmin\tavgLen\tavgStart\tavgOA\tavgQ\n";
	printf STAT ("%d\t%d\t%d\t%d\t%d\t%.4f\t%.0f\n",
		$remainQ[0],$max_bp[0],$min_bp[0],$sum_bp[0]/$remainQ[0],$sum_s[0]/$remainQ[0],$sum_oa[0]/$remainQ[0],$sum_Q[0]/$remainQ[0]);
	print STAT "Read_2\tmax\tmin\tavgLen\tavgStart\tavgOA\tavgQ\n";
	printf STAT ("%d\t%d\t%d\t%d\t%d\t%.4f\t%.0f\n",
		$remainQ[1],$max_bp[1],$min_bp[1],$sum_bp[1]/$remainQ[1],$sum_s[1]/$remainQ[1],$sum_oa[1]/$remainQ[1],$sum_Q[1]/$remainQ[1]);
	print STAT "Single\tmax\tmin\tavgLen\tavgStart\tavgOA\tavgQ\n";
	printf STAT ("%d\t%d\t%d\t%d\t%d\t%.4f\t%.0f\n",
		$remainQ[2],$max_bp[2],$min_bp[2],$sum_bp[2]/$remainQ[2],$sum_s[2]/$remainQ[2],$sum_oa[2]/$remainQ[2],$sum_Q[2]/$remainQ[2]);
}
print STAT "length\tcount\n";
foreach my $l(sort {$a<=>$b} keys %STAT){
	print STAT "$l\t$STAT{$l}\n";
}
close STAT;
exit;
# sub
sub fqRead {
	my $present= shift;
	my $handle = shift;
	chomp($_[0] = $present);
	chomp($_[1] = <$handle>);
	chomp($_[2] = <$handle>);
	chomp($_[3] = <$handle>);
	return(@_);
}

sub ca1_cut {
    my $Sc  = shift;
    my $cut = shift;
    my $sysQ = shift;
	my $seq = shift;
    my $q = shift;
    my ($ca, $min, $tmp,$oa0, $oa1, $s,$p) = (1, 0, 1, 1, 1, 0, 0);
    my @Q;
	my @PQ;
    # cal phrd Q
    while($p<length($q)){
        $_ = substr($q,$p,1);
        $_ = ord($_) - $sysQ;
		push @PQ, $_;
        $_ = 1 - 10**(-$_/10);
        $_ = ($_==0)?0.1:$_;
        push @Q, $_;
        $p ++;
    }
    # cal first seed OA
    $p = 0;
    my @seedOA = (1);
    while($p<30){
        $seedOA[0] *= $Q[$p];
        $p ++;
    }
    # choose best seed
    while($p<length($q)){
        $seedOA[$p-29] = $seedOA[$p-30] * $Q[$p] / $Q[$p-30];
        $s = ($seedOA[$p-29] > $seedOA[$s])?$p-29:$s;
        $p++;
        last if $seedOA[$s]  >= $Sc;
    }
    # trim
#    $ca = $seedOA[$s];
    $p = $s + 30;
    $ca = $seedOA[$s];
    while($p<length($q)){
        my $acc = $Q[$p];
#        $oa0 *= $acc;
        if($acc < $min){
            $oa1 = $ca * $min;
            $min = $acc;
        }else{
            $oa1 = $ca * $acc;
        }
        last if $oa1 < $cut;
        $p ++;
        $ca = $oa1;

    }
	my $PQ =0;
	for(my $i=$s;$i<$p;$i++){
		$PQ += $PQ[$i];
	}
	$p   = $p - $s;
	$seq = substr($seq,$s,$p);
	$q   = substr($q,$s,$p);
    return($ca,$PQ,$s,$p,$seq,$q);
}

sub cumulate {
	my $n = shift;
	my $m = shift;
	$remainQ[$n]  ++;
	$STAT{$len[$m]} ++;
	$max_bp[$n] = ($max_bp[$n] > $len[$m])?$max_bp[$n]:$len[$m];
	$min_bp[$n] = ($min_bp[$n] < $len[$m])?$min_bp[$n]:$len[$m];
	$sum_bp[$n] += $len[$m];
	$sum_s[$n]  += $start[$m];
	$sum_oa[$n] += $Aqual[$m];
	$sum_Q[$n] += $PQ[$m];
}



