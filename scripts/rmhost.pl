#!/usr/bin/perl -w

############################################################
#	Copyright (C) 2012 BGI-shenzhen-MicrobeGroup
#	Written by Yuanlin Guan (guanyuanlin@genomics.org.cn)
#	updated by fangchao@genomics.cn
#	Version: 1.2 (20160418)
############################################################
use strict;
use Getopt::Std;
use FindBin;
#use PerlIO::gzip; This pack frequently goes wrong.
sub openM {	#a sub of open file handle method
	my $f =shift;my $m =($f =~ /\.gz$/)?"gzip -dc $f|":"$f";
	return($m);
}

sub usage{
	print STDERR "usage: $0 <option> <value>...\n";
	print STDERR "Options\n";
	print STDERR "\t-a : read1 file, required\n";
	print STDERR "\t-b : read2 file, optional\n";
	print STDERR "\t-c : single read file, optional\n";
	print STDERR "\t-d : host database, required\n";
	print STDERR "\t-D : match mode, default 4\n";
	print STDERR "\t-s : seed length, default 30\n";
	print STDERR "\t-r : repeat hit, default 1\n";
	print STDERR "\t-m : minimal insert size, default 400\n";
	print STDERR "\t-x : maximal insert size, default 600\n";
	print STDERR "\t-n : filter low-quality reads containing >n Ns before alignment, [5]\n";
	print STDERR "\t-v : maximum number of mismatches, default 7\n";
	print STDERR "\t-i : identity, default 0.9\n";
	print STDERR "\t-t : number of processors, default 3\n";
	print STDERR "\t-f : simple soap result, default Y (Y/N)\n";
	print STDERR "\t-q : fastq reads\n";
	print STDERR "\t-p : output prefix, required\n";
	print STDERR "\t-h : show the help message\n";
	exit(1);
}

our ($opt_a, $opt_b, $opt_c, $opt_d, $opt_D, $opt_s, $opt_r, $opt_m, $opt_n);
our ($opt_x, $opt_v, $opt_i, $opt_t, $opt_f, $opt_q, $opt_p, $opt_h);
getopts('a:b:c:d:D:s:r:m:n:x:v:i:t:f:q:p:h');

&usage if $opt_h;
&usage unless $opt_a && $opt_d && $opt_p;
$opt_D = 4 unless $opt_D;
$opt_s = 30 unless $opt_s;
$opt_r = 1 unless $opt_r;
$opt_m = 400 unless $opt_m;
$opt_x = 600 unless $opt_x;
$opt_n = 5 unless $opt_n;
$opt_v = 7 unless $opt_v;
$opt_i = 0.9 unless $opt_i;
$opt_t = 3 unless $opt_t;
$opt_f = "Y" unless $opt_f;

my ($clean1, $clean2, $single, $stat);
my ($total_num, $total_single, $pair_num, $single_num, $read1_base, $read2_base, $single_base) = (0, 0, 0, 0, 0, 0, 0);
my ($read1_length, $read2_length, $max_length, $min_length, $avg_length);
$max_length = 0; $min_length = 1e10;

my ($soap_path, $shell, %remove, @host, $host);
my ($name1, $seq1, $head1, $qual1, $name2, $seq2, $head2, $qual2);

$soap_path = "$FindBin::Bin/soap2.22";
@host = split /:/, $opt_d;
$host = join " -D ", @host;

if($opt_b){
	if($opt_q){
		$clean1 = "$opt_p.rmhost.1.fq.gz";
		$clean2 = "$opt_p.rmhost.2.fq.gz";
		$single = "$opt_p.rmhost.single.fq.gz";
	}
	else{
		$clean1 = "$opt_p.rmhost.1.fa.gz";
		$clean2 = "$opt_p.rmhost.2.fa.gz";
		$single = "$opt_p.rmhost.single.fa.gz";
	}
}
else{
	if($opt_q){
		$clean1 = "$opt_p.rmhost.fq.gz";
	}
	else{
		$clean1 = "$opt_p.rmhost.fa.gz";
	}
}
$stat = "$opt_p.rmhost.stat_out";

if($opt_b){
	$shell = "$soap_path -a $opt_a -b $opt_b -D $host -M $opt_D -l $opt_s -r $opt_r -m $opt_m -x $opt_x -n $opt_n -v $opt_v -c $opt_i -p $opt_t ";
	if($opt_f eq "Y"){
		$shell .= "-S ";
	}
	$shell .= "-o $opt_p.rmhost.soap.pe -2 $opt_p.rmhost.soap.se 2> $opt_p.rmhost.soap.log";
	if(system($shell)){
		print STDERR "reads align host error\n";
		exit(1);
	}
	&get_remove("$opt_p.rmhost.soap.pe");
	&get_remove("$opt_p.rmhost.soap.se");
	if($opt_q){
		&get_paired_fq($opt_a, $opt_b, $clean1, $clean2);
	}
	else{
		&get_paired_fa($opt_a, $opt_b, $clean1, $clean2);
	}
	if($opt_c){
		$shell = "$soap_path -a $opt_c -D $host -M $opt_D -l $opt_s -r $opt_r -n $opt_n -v $opt_v -c $opt_i -p $opt_t ";
		if($opt_f eq "Y"){
			$shell .= "-S ";
		}
		$shell .= "-o $opt_p.rmhost.soap.single 2> $opt_p.rmhost.soap.single.log";
		if(system($shell)){
			print STDERR "single read align host error\n";
			exit(1);
		}
		&get_remove("$opt_p.rmhost.soap.single");
		if($opt_q){
			&get_single_fq($opt_c, $single);
		}
		else{
			&get_single_fa($opt_c, $single);
		}
	}
}
else{
    $soap_path='soap2.2d' if -x 'soap2.2d';
	$shell = "$soap_path -a $opt_a -D $host -M $opt_D -l $opt_s -r $opt_r -n $opt_n -v $opt_v -c $opt_i -p $opt_t ";
	if($opt_f eq "Y"){
		$shell .= "-S ";
	}
    if(-x 'soap2.2d'){
        $shell .= "-o $opt_p.rmhost.soap -u $clean1 2> $opt_p.rmhost.soap.log";
        if(system($shell)){
            print STDERR "single read align host error\n";
            exit(1);
        }
		system("gzip $clean1");
	
    } else {
    	$shell .= "-o $opt_p.rmhost.soap 2> $opt_p.rmhost.soap.log";
    	if(system($shell)){
    		print STDERR "single read align host error\n";
    		exit(1);
    	}
    	&get_remove("$opt_p.rmhost.soap");
    	if($opt_q){
    		&get_single_fq($opt_a, $clean1);
    	}
    	else{
    		&get_single_fa($opt_a, $clean1);
    	}
    }
}

open O, ">$stat" or die "can't open file $stat $!\n";
if($opt_b){
	$avg_length = ($read1_base + $read2_base + $single_base) / ($pair_num * 2 + $single_num);
	print O "Total_pair\tAligned_pair\tTotal_single\tAligned_single\tread1_base\tread2_base\tsingle_base\tmax_length\tmin_length\tavg_length\n";
	print O "$total_num\t$pair_num\t$total_single\t$single_num\t$read1_base\t$read2_base\t$single_base\t$max_length\t$min_length\t$avg_length\n";
}
else{
	$avg_length = $single_base / $single_num;
	print O "Total_single\tAligned_single\tsingle_base\tmax_length\tmin_length\tavg_length\n";
	print O "$total_single\t$single_num\t$single_base\t$max_length\t$min_length\t$avg_length\n";
}
close O;
#system("rm -f $opt_p.rmhost.soap $opt_p.rmhost.soap.single $opt_p.rmhost.soap.pe $opt_p.rmhost.soap.se");


sub get_remove{
	my $file = shift;
	open I, &openM($file) or die "can't open file $file $!\n";
	while(<I>){
		chomp;
		my @temp = split /\s+/, $_;
		$temp[0] =~ s/\/[12]$//;
		$remove{$temp[0]} = undef;
	}
}

sub get_paired_fq{
	my ($in1, $in2, $out1, $out2) = @_;
	open I1, &openM($in1) or die "can't open file $in1 $!\n";
	open I2, &openM($in2) or die "can't open file $in2 $!\n";
	open O1, "|gzip >$out1" or die "can't open file $out1 $!\n";
	open O2, "|gzip >$out2" or die "can't open file $out2 $!\n";

	while(<I1>){
		$total_num++;
		chomp;
		my @heads = split /\s+/, $_;
        $name1 = $heads[0];
		$name1 =~ s/^\@//;
		$name1 =~ s/\/[12]//;
		
		$name2 = <I2>;
		chomp $name2;
        my @heads2 = split /\s+/, $name2;
        $name2 = $heads2[0];
		$name2 =~ s/^\@//;
		$name2 =~ s/\/[12]//;
		die "reads file error. Line $.: name1:$name1 ne name2:$name2\n" if $name1 ne $name2;
		if(exists $remove{$name1}){
			<I1>;<I1>;<I1>;
			<I2>;<I2>;<I2>;
		}
		else{
			$seq1 = <I1>; $head1 = <I1>; $qual1 =<I1>;
			$seq2 = <I2>; $head2 = <I2>; $qual2 = <I2>;
			print O1 "\@$name1/1\n$seq1$head1$qual1";
			print O2 "\@$name2/2\n$seq2$head2$qual2";
			$pair_num++;
			chomp $seq1; chomp $seq2;
			$read1_length = length($seq1);
			$read2_length = length($seq2);
			$read1_base += $read1_length;
			$read2_base += $read2_length;
			$max_length = $max_length > $read1_length ? $max_length : $read1_length;
			$min_length = $min_length < $read1_length ? $min_length : $read1_length;
			$max_length = $max_length > $read2_length ? $max_length : $read2_length;
			$min_length = $min_length < $read2_length ? $min_length : $read2_length;
		}
	}
	close I1;
	close I2;
	close O1;
	close O2;
}

sub get_paired_fa{
	my ($in1, $in2, $out1, $out2) = @_;
	open I1, &openM($in1) or die "can't open file $in1 $!\n";
	open I2, &openM($in2) or die "can't open file $in2 $!\n";
	open O1, "|gzip >", "$out1" or die "can't open file $out1 $!\n";
	open O2, "|gzip >", "$out2" or die "can't open file $out2 $!\n";

	$/ = ">";
	<I1>; <I2>;
	while(<I1>){
		$total_num++;
		chomp;
		my @temp = split /\n/, $_, 2;
		$name1 = $temp[0];
		$name1 =~ s/\/[12]$//;
		$seq1 = $temp[1];
		$seq1 =~ s/\n//g;
		my $tmp = <I2>;
		chomp $tmp;
		@temp = split /\n/, $_, 2;
		$name2 = $temp[0];
		$name2 =~ s/\/[12]$//;
		$seq2 = $temp[1];
		$seq2 =~ s/\n//g;
		die "reads file error. Line $.: name1:$name1 ne name2:$name2\n" if $name1 ne $name2;
		unless(exists $remove{$name1}){
			print O1 "\>$name1/1\n$seq1\n";
			print O2 "\>$name2/2\n$seq2\n";
			$pair_num++;
			$read1_length = length($seq1);
			$read2_length = length($seq2);
			$read1_base += $read1_length;
			$read2_base += $read2_length;
			$max_length = $max_length > $read1_length ? $max_length : $read1_length;
			$min_length = $min_length < $read1_length ? $min_length : $read1_length;
			$max_length = $max_length > $read2_length ? $max_length : $read2_length;
			$min_length = $min_length < $read2_length ? $min_length : $read2_length;
		}
	}
	$/ = "\n";
	close I1;
	close I2;
	close O1;
	close O2;
}

sub get_single_fq{
	my ($in, $out) = @_;
	open I, &openM($in) or die "can't open file $in $!\n";
	open O, "|gzip > $out" or die "can't open file $out $!\n";

	while(<I>){
		$total_single++;
		chomp;
		my @temp = split /\s+/, $_;
		$name1 = $temp[0];
		my $name = $temp[0];
		$name1 =~ s/^\@//;
		$name1 =~ s/\/[12]$//;
		if(exists $remove{$name1}){
			<I>;<I>;<I>;
		}
		else{
			$seq1 = <I>; $head1 = <I>; $qual1 =<I>;
			print O "$name\n$seq1$head1$qual1";
			$single_num++;
			chomp $seq1;
			$read1_length = length($seq1);
			$single_base += $read1_length;
			$max_length = $max_length > $read1_length ? $max_length : $read1_length;
			$min_length = $min_length < $read1_length ? $min_length : $read1_length;
		}
	}
	close I;
	close O;
}

sub get_single_fa{
	my ($in, $out) = @_;
	open I, &openM($in) or die "can't open file $in $!\n";
	open O, "|gzip > $out" or die "can't open file $out $!\n";

	$/ = ">";
	<I>;
	while(<I>){
		$total_single++;
		chomp;
		my @temp = split /\s+/, $_, 2;
		$name1 = $temp[0];
		my $name = $temp[0];
		$name1 =~ s/\/[12]$//;
		$seq1 = $temp[1];
		$seq1 =~ s/\n//g;
		unless(exists $remove{$name1}){
			print O "\>$name\n$seq1\n";
			$single_num++;
			$read1_length = length($seq1);
			$single_base += $read1_length;
			$max_length = $max_length > $read1_length ? $max_length : $read1_length;
			$min_length = $min_length < $read1_length ? $min_length : $read1_length;
		}
	}
	close I;
	close O;
}

