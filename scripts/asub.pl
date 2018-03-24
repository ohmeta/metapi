#!/usr/bin/perl -w

# Author: lh3

use strict;
use warnings;
use Getopt::Std;
use Cwd qw/abs_path/;
use Sys::Hostname qw/hostname/;

my $version = "2.1";

# parsing command line
my %opts = (R=>'', j=>'', q=>'', c=>'', w=>'', M=>'', G=>'', n=>1, m=>'', W=>'', C=>'', g=>1, P=>'smp');
getopts('P:R:j:k:q:c:pw:M:n:m:W:g:C:G', \%opts);
&usage if (-t STDIN && @ARGV == 0);

# adjusting parameters
$opts{g} = 1 if $opts{g} < 1;
$opts{n} = $opts{g} if !$opts{G} && $opts{n} < $opts{g};
$opts{j} ||= sprintf("asub_$$%.3d", rand(1000));

if (!defined($opts{k}) || $opts{k} < 1) { # prepare to run
	# LSF or Grid Engine?
	my $is_slurm = &which('sbatch')? 1 : 0;
	my $is_ge = &which('bsub')? 0 : 1;
	die ("ERROR: option -g only works with LSF for now\n") if ($is_ge && $opts{g} != 1);

	my $n_cmds = 0;
	open(FILE, ">$opts{j}.sh") || die;
	while (<>) {
		s/\s*&\s*$//; # remove "&" at the end of the command
		next if /^\s*$/;
		print FILE $_;
		++$n_cmds;
	}
	close(FILE);
	my $n_jobs = int(($n_cmds + $opts{g} - 1) / $opts{g});

	my $cmd;
	$opts{c} = $n_jobs if $opts{c} && $opts{c} > $n_jobs;
	$opts{G} = $opts{G}? "-G" : "";
	$_ = $opts{j};
	if ($is_slurm) {
		$opts{q} = qq/-p $opts{q}/ if $opts{q};
		$opts{R} = qq/--mem=/ . ($opts{M}*1024) if $opts{M};
		$opts{R} .= qq/ -N 1 -n $opts{n}/ if $opts{n} > 1;
		$opts{R} .= qq/ -t $opts{W}/ if $opts{W};
		$cmd = qq(mkdir -p $_.out $_.err && echo -e '#!/bin/sh\\n$0 $opts{G} -g $opts{g} -k \${SLURM_ARRAY_TASK_ID} $_.sh' | sbatch -J$_ )
			. qq(--array=1-$n_jobs -o $_.out/\%a.out -e $_.err/\%a.err $opts{R} $opts{q});
	} elsif (!$is_ge) { # LSF
		$opts{q} = qq/-q $opts{q}/ if $opts{q};
		$opts{c} = qq/\%$opts{c}/ if $opts{c};
		$opts{w} = qq/-w "$opts{w}"/ if $opts{w};

		# resource requirements and limits
		$opts{R} = qq/-R "$opts{R}"/ if $opts{R};
		$opts{R} .= qq/ -m $opts{m}/ if $opts{m};
		$opts{R} .= qq/ -n $opts{n}/ if $opts{n};
		$opts{R} .= qq/ -W $opts{W}/ if $opts{W};
		$opts{R} .= qq/ -c $opts{C}/ if $opts{C};
		$opts{R} .= qq/ -R "span[hosts=1]"/ if $opts{n};
		if ($opts{M}) {
			my $mem = int($opts{M} * 1024 / $opts{n}); # in LSF, rusage and -M set per-process limit
			$opts{R} .= qq/ -R "rusage[mem=$mem]"/;
			#$opts{R} .= qq/ -M $mem/;
		}

		# the command line to bsub
		$cmd = qq(mkdir -p $_.out $_.err && echo '$0 $opts{G} -g $opts{g} -k \${LSB_JOBINDEX} $_.sh' | bsub -J$_)
			   . qq("[1-$n_jobs]$opts{c}" -o $_.out/\%I.out -e $_.err/\%I.err $opts{R} $opts{q} $opts{w});
	} else { # SGE
		$opts{q} = qq/-q $opts{q}/ if $opts{q};
		$opts{R} = qq/-l "$opts{R}"/ if $opts{R};
		$opts{R} .= qq/ -pe $opts{P} $opts{n}/ if $opts{n} > 1;
		$opts{R} .= qq/ -l h_rt=$opts{W} -l s_rt=$opts{W}/ if $opts{W};
		if ($opts{M}) {
			my $mem = int($opts{M} * 1024 / $opts{n});
			$opts{R} .= qq/ -l h_vmem=$mem/ . qq/M/;
		}
		$opts{w} = qq/ -hold_jid $opts{w}/ if $opts{w};
		my $asub = abs_path($0);
		$cmd = qq(mkdir -p $_.out $_.err && echo '$asub $opts{G} -g $opts{g} -k \${SGE_TASK_ID} $_.sh' | qsub -V -N $_ -cwd )
			. qq(-t 1-$n_jobs -o '$_.out/\$TASK_ID.out' -e '$_.err/\$TASK_ID.err' $opts{R} $opts{q} $opts{w});
	}
	defined($opts{p})? (print "$cmd\n") : system($cmd);
} else { # run the command
	warn("[asub] Hostname: ".hostname.", Arch-OS: ".&get_cpu_sys."\n");
	my $lineno = 0;
	my @cmds = ();
	while (<>) {
		chomp;
		++$lineno;
		push(@cmds, $_) if $lineno > ($opts{k} - 1) * $opts{g} && $lineno <= $opts{k} * $opts{g};
	}

	# launch command lines
	my ($pid, %pids);
    my $shell = $ENV{'SHELL'} || '/bin/sh';
	warn("[asub] --- BEGIN OF COMMAND STDERR ---\n");
	if ($opts{G}) { # in serial
		$pid = 0;
		for my $cmd (@cmds) {
			++$pid;
			$pids{"$pid"}[0] = $cmd;
			system("$shell", "-c", $cmd);
			$pids{"$pid"}[1] = $? >> 8;
		}
	} else { # in parallel
		for my $cmd (@cmds) {
			if (!defined($pid = fork())) { # fork failure
				die("[asub] failed to fork. Abort!\n");
			} elsif ($pid != 0) { # this the parent process
				$pids{$pid}[0] = $cmd;
			} else { # this is the child process
				exec "$shell", "-c", $cmd;
			}
		}
		while (($pid = wait()) >= 0) {
			$pids{$pid}[1] = $? >> 8;
		}
	}
	warn("[asub] --- END OF COMMAND STDERR ---\n");

	# print return value for each command line
	my $ret = 0;
	for my $pid (sort {$a<=>$b} keys %pids) {
		$ret = $pids{$pid}[1] if $pids{$pid}[1] != 0;
		warn("[asub] The following command (PID=$pid) returns $pids{$pid}[1]\n");
		warn("[asub]   $pids{$pid}[0]\n");
	}
	exit $ret;
}

sub get_cpu_sys {
	my $dir = `uname -p 2>/dev/null`;
	$dir = `uname -m 2>/dev/null` if (!$dir || $dir =~ "unknown");
	$dir .= '-'.`uname -s`;
	$dir = lc($dir);
	$dir =~ s/\s//g;
	return $dir;
}

sub which {
	my $file = shift;
	my $path = (@_)? shift : $ENV{PATH};
	return if (!defined($path));
	foreach my $x (split(":", $path)) {
		$x =~ s/\/$//;
		return "$x/$file" if (-x "$x/$file");
	}
	return;
}

sub usage {
  die qq(
Program: asub (Array Job submission with bsub/qsub)
Version: $version
Contact: Heng Li <lh3\@sanger.ac.uk>\n
Usage:   asub [options] <cmd_file>\n
Options: -R STR    resources string (only one -R allowed) [null]
         -M INT    maximum total memory in GB [null]
         -W STR    runtime limit [none]
         -n INT    #CPUs per job [$opts{n}]
         -C STR    CPU time limit (LSF only) [none]
         -w STR    job dependency [null]
         -q STR    queue name [default queue]
         -j STR    job name [auto]
         -c INT    number of concurrent jobs (LSF only) [max]
         -m STR    host group (LSF only) [null]
         -g INT    group INT command lines in one job [$opts{g}]
         -G        serializing {-g}
         -p        print the submission command only

Note: For option -R and -w, SGE and LSF have different syntax. Please
      check SGE/LSF manual page for details. Here are some examples for
      LSF:

        -R "select[type==X86_64&&mem>800] rusage[mem=800]"
        -w "done(my_job_001)"

      And for SGE:

        -R h_cpu=86400,h_data=1000000000
        -w my_job_001

);
}
