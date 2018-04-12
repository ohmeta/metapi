#!/hwfssz1/ST_META/CD/zhujie/program/bioenv/bin/perl -w
# =============================================================================
# CD-HIT
# http://cd-hit.org/
# http://bioinformatics.burnham.org/cd-hit
#
# program written by
#                                      Weizhong Li
#                                      UCSD, San Diego Supercomputer Center
#                                      La Jolla, CA, 92093
#                                      Email liwz@sdsc.edu
# modified by Weineng Chen
# =============================================================================

use strict;
use Cwd;
no strict "refs";

my $script_name = $0; # script cd-hit-para.v2.pl absolute path
my $script_dir = $0;
$script_dir =~ s/[^\/]+$//; # remove last / symbol
$script_dir = "./" unless ($script_dir); # 如果没有$script_dir，则把$script_dir设置为当前目录
chop($script_dir); # 去掉$script_dir末尾的空格

my $arg; # 参数
my $in;  # 输入的fasta，要去冗余
my $out; # 输出目录
my $arg_pass          = "";
my $para              = "";
my $host_no           = 0;
my @hosts             = ();

# 一些需要调用的程序的absolute path
my $cd_hit_div_exe    = "$script_dir/cd-hit-div";
my $cd_hit_exe        = "$script_dir/cd-hit";
my $cd_hit_est_exe    = "$script_dir/cd-hit-est";
my $cd_hit_2d_exe     = "$script_dir/cd-hit-2d";
my $cd_hit_est_2d_exe = "$script_dir/cd-hit-est-2d";
my $clstr_merge_exe   = "$script_dir/clstr_merge.pl";

my $seg_no            = 64;
my $restart_in        = "";

my $queue             = 0;
my $local_cpu         = 0;
my $queue_type        = "PBS";
my $prog              = "cd-hit";

# 这个stable_files是干什么用的呢？
my %stable_files      = ();


# 投递任务到SGE集群需要提供的信息
# 队列名称，项目ID，向SGE资源调度系统申请的资源(内存和CPU数目)
my $queue_list        = "";
my $project_name      = "";
my $resource          = "";

# cd-hit-para.v2.pl
#    -i 292S_oral_allgene.fa # 输入的fasta
#    -o 292S_oral_nr # 输出的前缀
#    --S 64 # 切分成64分子fasta
#    --Q 90 # 最多投递90个任务
#    --T SGE # SGE投递平台
#    -q st.q # 队列名称
#    -P F16ZQSB1SY2779 # 项目ID
#    --l vf=4G,p=4 # 申请的资源
#    --R 292S_oral_nr.restart # restart文件 # 第一次跑的时候不要指定
#    --P cd-hit-est # 调用的实际的去冗余程序

#    # 接下来的参数都是cd-hit-est的参数
#    -G 0 # 采用局部比对
#    -n 8 # word size设置为8
#    -aS 0.9 # 短reads长度要占到长reads的90%
#    -c 0.95 # 序列比对的identity
#    -M 13000 # 内存限制，最多用13G
#    -d 0 
#    -r 1 
#    -g 1 
#    -T 4 # 4个线程
=pod
====== CD-HIT version 4.6 (built on Sep 18 2016) ======

Usage: cd-hit-est [Options]

Options

   -i   input filename in fasta format, required
   -o   output filename, required
   -c   sequence identity threshold, default 0.9
        this is the default cd-hit's "global sequence identity" calculated as:
        number of identical amino acids in alignment
        divided by the full length of the shorter sequence
   -G   use global sequence identity, default 1
        if set to 0, then use local sequence identity, calculated as :
        number of identical amino acids in alignment
        divided by the length of the alignment
        NOTE!!! don't use -G 0 unless you use alignment coverage controls
        see options -aL, -AL, -aS, -AS
   -b   band_width of alignment, default 20
   -M   memory limit (in MB) for the program, default 800; 0 for unlimitted;
   -T   number of threads, default 1; with 0, all CPUs will be used
   -n   word_length, default 10, see user's guide for choosing it
   -l   length of throw_away_sequences, default 10
   -d   length of description in .clstr file, default 20
        if set to 0, it takes the fasta defline and stops at first space
   -s   length difference cutoff, default 0.0
        if set to 0.9, the shorter sequences need to be
        at least 90% length of the representative of the cluster
   -S   length difference cutoff in amino acid, default 999999
        if set to 60, the length difference between the shorter sequences
        and the representative of the cluster can not be bigger than 60
   -aL  alignment coverage for the longer sequence, default 0.0
        if set to 0.9, the alignment must covers 90% of the sequence
   -AL  alignment coverage control for the longer sequence, default 99999999
        if set to 60, and the length of the sequence is 400,
        then the alignment must be >= 340 (400-60) residues
   -aS  alignment coverage for the shorter sequence, default 0.0
        if set to 0.9, the alignment must covers 90% of the sequence
   -AS  alignment coverage control for the shorter sequence, default 99999999
        if set to 60, and the length of the sequence is 400,
        then the alignment must be >= 340 (400-60) residues
   -A   minimal alignment coverage control for the both sequences, default 0
        alignment must cover >= this value for both sequences
   -uL  maximum unmatched percentage for the longer sequence, default 1.0
        if set to 0.1, the unmatched region (excluding leading and tailing gaps)
        must not be more than 10% of the sequence
   -uS  maximum unmatched percentage for the shorter sequence, default 1.0
        if set to 0.1, the unmatched region (excluding leading and tailing gaps)
        must not be more than 10% of the sequence
   -U   maximum unmatched length, default 99999999
        if set to 10, the unmatched region (excluding leading and tailing gaps)
        must not be more than 10 bases
   -B   1 or 0, default 0, by default, sequences are stored in RAM
        if set to 1, sequence are stored on hard drive
        it is recommended to use -B 1 for huge databases
   -p   1 or 0, default 0
        if set to 1, print alignment overlap in .clstr file
   -g   1 or 0, default 0
        by cd-hit's default algorithm, a sequence is clustered to the first
        cluster that meet the threshold (fast cluster). If set to 1, the program
        will cluster it into the most similar cluster that meet the threshold
        (accurate but slow mode)
        but either 1 or 0 won't change the representatives of final clusters
   -r   1 or 0, default 1, by default do both +/+ & +/- alignments
        if set to 0, only +/+ strand alignment
   -mask        masking letters (e.g. -mask NX, to mask out both 'N' and 'X')
   -match       matching score, default 2 (1 for T-U and N-N)
   -mismatch    mismatching score, default -2
   -gap gap opening score, default -6
   -gap-ext     gap extension score, default -1
   -bak write backup cluster file (1 or 0, default 0)
   -h   print this help

   Questions, bugs, contact Limin Fu at l2fu@ucsd.edu, or Weizhong Li at liwz@sdsc.edu
   For updated versions and information, please visit: http://cd-hit.org

   cd-hit web server is also available from http://cd-hit.org

   If you find cd-hit useful, please kindly cite:

   "Clustering of highly homologous sequences to reduce thesize of large protein database", Weizhong Li, Lukasz Jaroszewski & Adam Godzik. Bioinformatics, (2001) 17:282-283
   "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik. Bioinformatics, (2006) 22:1658-1659
=cut

# 解析参数
while ($arg = shift) 
{
    if    ($arg eq "-h" ) { print_usage(); }
    elsif ($arg eq "-i" ) { $in         = shift; }
    elsif ($arg eq "-o" ) { $out        = shift; }
    elsif ($arg eq "--B") { $para       = shift; }
    elsif ($arg eq "--L") { $local_cpu  = shift; }
    elsif ($arg eq "--P") { $prog       = shift; }
    elsif ($arg eq "--S") { $seg_no     = shift; }
    elsif ($arg eq "--Q") { $queue      = shift; }
    elsif ($arg eq "--T") { $queue_type = shift; }
    elsif ($arg eq "--R") { $restart_in = shift; }
    elsif ($arg eq "-q")  { $queue_list = shift; }
    elsif ($arg eq "-P")  { $project_name = shift; }
    elsif ($arg eq "--l") { $resource = shift; }
    # 不在 -h -i -o --B --L --P --S --Q --T --R -q -P --l的参数都连接到$arg_pass这个字符串变量里面
    else  { $arg_pass .= " $arg " . shift; }
}
# 没有指定输入和输出前缀就打印程序usage
($in and $out) || print_usage();

# 设置具体的调用程序
# 如果是处理DNA，就把cd_hit_ext指向cd_hit_est_exe
# 把cd_hit_2d_exe指向cd_hit_est_2d_exe
if ($prog eq "cd-hit-est") 
{
    $cd_hit_exe        = $cd_hit_est_exe;
    $cd_hit_2d_exe     = $cd_hit_est_2d_exe;
}

#my $pwd = `pwd`; chomp($pwd);

# 获取执行这个cd-hit-para.v2.pl时用户所在的目录
my $pwd = getcwd;
#die "$pwd\n";

# 设置工作目录，从这里看出$out只是输出前缀
# 但是最终的确会生成一个$out的文件，在$pwd目录下面
# 所以$outfile是reduced之后的结果
my $work_dir       = "$out.cd-hit-para-tmp";
my $infile         = (split /\//, $in)[-1];  # basename 292S_oral_allgene.fa
my $outfile        = (split /\//, $out)[-1]; # basename 292S_oral_nr

# restart文件，restart是这个脚本的精华
my $restart_file   = "$out.restart";          # 292S_oral_nr.restart
my $indiv          = "$work_dir/$infile.div"; # 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div

# 执行命令都放在一个数组里面
my @commands       = ();
# 每个命令的状态
my @command_status = ();
my $command_no     = 0;
my $cmd;
my ($i, $j, $k, $i1, $j1, $k1);

# 文件锁目录
my $lock_dir       = "$work_dir/$outfile-lock";      # 292S_oral_nr.cd-hit-para-tmp/292S_oral_nr-lock
# qsub日志目录
my $qsub_log       = "$work_dir/$outfile-qsub_log";  # 292S_oral_nr.cd-hit-para-tmp/292S_oral_nr-qsub_log


# readin a list of hosts
# 指定一系列的主机名称，然后添加到hosts数组里面
# 我们用SGE集群，一般不指定具体的主机
if ($para) 
{
    open(PARA, "$para") || die "can not open $para";
    while(my $ll = <PARA>) 
    {
        chop($ll); 
        $ll =~ s/\s//g;
        next unless ($ll);
        push(@hosts, $ll); 
        $host_no++; # 对主机数量进行计数
    }
    close(PARA);
}

# 比如 $queue = 90
if ($queue) 
{
    for ($i = 0; $i < $queue; $i++) 
    {
        push(@hosts, "queue_host.$i");
    }
    $host_no = $queue; # 如果指定了--Q, $host_no被赋值为$queue
}

# 如果程序在本地跑，设置其要用的CPU数目
if ($local_cpu) 
{
    for ($i = 0; $i < $local_cpu; $i++) 
    {
        push(@hosts, "localhost.$i");
    }
    $host_no = $local_cpu;
}
die "no host" unless $host_no;

# 如果指定了$restart_in文件，就 调用read_restart()这个函数
if (-e $restart_in) 
{
    read_restart();
}
# 不然, 如果相应的目录不存在，就创建相应的目录
else 
{
    $cmd = `mkdir $work_dir` unless (-e $work_dir); # mkdir 292S_oral_nr.cd-hit-para-tmp
    `mkdir $lock_dir` unless (-e $lock_dir);        # mkdir 292S_oral_nr.cd-hit-para-tmp/292S_oral_nr-lock
    `mkdir $qsub_log` unless (-e $qsub_log);        # mkdir 292S_oral_nr.cd-hit-para-tmp/292S_oral_nr-qsub_log
    # 分配执行命令
    assign_commands();
    # 刷新restart文件
    write_restart();
}

#exit(1);
#dbdiv run on master node?
if ($command_status[0] eq "wait") 
{
    # 任务状态wait，执行该任务
    $cmd = `$commands[0]`;
    # 执行完后设置任务状态为done
    $command_status[0] = "done";
    # 刷新restart文件
    write_restart();

    # 这步是在做什么呢？
    for ($i=0; $i<$seg_no; $i++) 
    {
        my $idb = "$indiv-$i";
        $stable_files{$idb} = 1;
    }
}
#main runing loop
my $sleep_time = 1;
while(1) 
{
#refresh job status by checking output files
#check whether all jobs are done or not
    my $finish_flag = 1;
    my $status_change = 0;

    for ($i = 1; $i < $command_no; $i++) 
    {
        next if ($command_status[$i] eq "done");
        $finish_flag = 0;
        next if ($command_status[$i] eq "wait");
        my $tcmd = $commands[$i];
        my $output = "";
        if ($tcmd =~ / -o\s+(\S+)/)
        {
            $output = $1;
            # 检查执行完后的任务的输出文件是否存在
            if ((-s $output) && (-s "$output.clstr") && (-s "$output.done")) 
            {
                # 如果存在则改变任务状态，已完成
                $command_status[$i] = "done";
                $status_change = 1;
            }
#change some error run to wait status, run again.
#            elsif (-s $output){
#                wait_stable_file($output);
#                wait_stable_file("$output.clstr");
#                unless (-s "$output.clstr"){
#                    $stable_files{"$output.clstr"} = 0;
#                    $command_status[$i] = "wait";
#                    $status_change = 1;
#                }
#            }
#            elsif(-s "$output.done"){
#                $command_status[$i] = "wait";
#                `rm -f $output.done`;
#                $status_change = 1;
#            }
        }
    }
    if ($status_change) 
    {
        write_restart();
    }
    else 
    {
        sleep($sleep_time); print ".";
    }
    last if $finish_flag;

    my $job_sent = 0;
    for ($i=1; $i<$command_no; $i++) 
    {
        # 已完成的任务不再投递
        next if ($command_status[$i] eq "done");
        # 已投递上正在运行的任务也不再投递
        next if ($command_status[$i] eq "run");
        my $tcmd = $commands[$i];
        my $in1 = "";
        my $in2 = "";
        my $out_done = "";
        if ($tcmd =~ / -i\s+(\S+)/ ) {$in1 = $1;}
        if ($tcmd =~ / -i2\s+(\S+)/) {$in2 = $1;}
        if ($tcmd =~ / -o\s+(\S+)/ ) {$out_done = "$1.done";}
        my $input_flag = 0;
#    die "$in1\t$in2\n"

        if (($in1 =~ /\S/) and ($in2 =~ /\S/)) 
        {
            $input_flag = 1 if ((-s $in1) and (-s "$in1.clstr") and (-s "$in1.done") and (-s $in2) and (-s "$in2.clstr") and (-s "$in2.done"));
            $input_flag = 1 if ((-s $in1) and (-s "$in1.clstr") and (-s "$in1.done") and (-s $in2) and ($in2 =~ /(div-\d+)$/));
        }
        elsif ($in1 =~ /\S/) 
        {
            $input_flag = 1 if ((-s $in1) and (-s "$in1.clstr") and (-s "$in1.done"));
            $input_flag = 1 if ((-s $in1) and ($in1 =~ /(div-\d+)$/));
        }
        else 
        {
            die "Error at $tcmd\n";
        }
        next unless $input_flag;

#now input files are ready, wait
        wait_stable_file($in1);
        wait_stable_file($in2) if ($in2 =~ /\S/);

        my $thost_idx = wait_for_available_host();
        my $thost     = $hosts[$thost_idx];
        my $tsh   = "$lock_dir/para-$thost_idx.sh";
        my $tlock = "$lock_dir/para-$thost_idx.lock";
        my $trm   = "";
##       $trm   = "rm -f $in2" if ($in2 =~ /\S/);
       do {$trm = "rm -f $in2"; $tcmd.=" && $trm"; } if ($in2 =~ /\S/);

        open(TSH, "> $tsh") || die;
        $cmd = `date > $tlock`;
        print TSH <<EOD;
#!/bin/sh
date > $tlock
$tcmd
##$trm
rm -f $tlock
date > $out_done
EOD
        close(TSH);
        if ($local_cpu) {
            $cmd = `sh $tsh >/dev/null 2>&1 &`;
            $command_status[$i] = "run";
            print "run at $thost $tcmd\n";
        }
        elsif ($queue and ($queue_type eq "PBS")) {
            my $t = "para-$thost_idx";
            open(QUEUE,"| qsub -N $t -o $qsub_log/$i.log -o $qsub_log/$i.log");
            print QUEUE "cd $pwd; sh $tsh";
            close(QUEUE);
            $command_status[$i] = "run";
        }
        # SGE 集群投递方式
        elsif ($queue and ($queue_type eq "SGE")) {
            $pwd =~ s/^\.//;
        #   die "$pwd\t$tsh\n";# rm 
            print STDERR "Submit job\n";	
            my $t = "para-$thost_idx";
            my $parameter="";
            $parameter.="-q $queue_list " if ($queue_list);
            $parameter.="-P $project_name " if ($project_name);
            $parameter.="-l $resource " if ($resource);
            open(QUEUE,"| qsub -cwd -N $t -e $qsub_log/$i.log -o $qsub_log/$i.log $parameter");
            print QUEUE <<EOD;
cd $pwd
sh $tsh

EOD
            close(QUEUE);
            print STDERR "Finish submitting job\n";

            $command_status[$i] = "run";
        }
        else {
            $cmd = `ssh -xqf $thost 'cd $pwd; sh $tsh  >/dev/null 2>&1 &'`;
            $command_status[$i] = "run";
            print "run at $thost $tcmd\n";
        }
        $sleep_time = 1;
        # 任务已sent
        $job_sent = 1;
        last;
    }

    if ((not $job_sent) and ($sleep_time < 60)) {
        $sleep_time +=5;
    }
    if ($job_sent) { write_restart();}
    print "*";

} ############ main run loop 

######## merge all .clstr file
my $out_clstr = "$out.clstr";
if (not -s $out_clstr) {

    my @reps = ();
    for ($i = 0; $i < $seg_no; $i++) {
        my $master_clstr = "$indiv-$i-/o.clstr";
        die "No file $master_clstr\n" unless (-s $master_clstr);

        my $this_rep = "$indiv-$i-/o";
        die "No rep $this_rep\n" unless (-e $this_rep);
        push(@reps, $this_rep);

        my @slave_clstr = ();
        for ($j = $i+1; $j < $seg_no; $j++) {
            my $tclstr = "$indiv-$j-/vs.$i.clstr";
            if (-s $tclstr) { push(@slave_clstr,$tclstr); }
            else {die "No file $tclstr\n";}
        }

        if (@slave_clstr) {
            my $tclstrs = join(" ", @slave_clstr);
            print  "$clstr_merge_exe $master_clstr $tclstrs >> $out_clstr\n";
            $cmd = `$clstr_merge_exe $master_clstr $tclstrs >> $out_clstr`;
        }
        else { #this is the last piece
            print  "cat $master_clstr >> $out_clstr";
            $cmd = `cat $master_clstr >> $out_clstr`;    
        }
    }

    my $out_clstr_ren = "$out.clstr.$$";
    open(TMP, $out_clstr) || die;
    open(OTMP, "> $out_clstr_ren") || die;
    my $no = 0;
    my $cno;
    my $ll;
    while($ll=<TMP>){
        if ($ll =~ /^>Cluster (\d+)/) {
            print OTMP ">Cluster $no\n"; $no++;
            $cno  = 0;
        }
        else {
            $ll =~ s/^\d+/$cno/;
            print OTMP $ll;
            $cno++;
        }
    }
    close(TMP);
    close(OTMP);
    sleep(10);
    $cmd = `mv $out_clstr_ren $out_clstr`;

    my $reps = join(" ", @reps);
    $cmd = `cat $reps > $out`;
}

if (1) {
    $cmd = `grep CPU $work_dir/*/log`;
    my @lls = split(/\n/, $cmd);
    my $cpu = 0;
    my $ll;
    foreach $ll (@lls) {
        if ($ll =~ /CPU\s+time\s+(\d+)/) {
            $cpu += $1;
        }
    }
    print "Total CPU time: $cpu\n";
}


sub wait_for_available_host {
    my ($i, $j, $k);
    my $sleep = 30;
    while(1) {

        for ($i=0; $i<$host_no; $i++) {
            my $thost = $hosts[$i];
            my $tlock = "$lock_dir/para-$i.lock";
            next if (-e $tlock);
            return $i;
        }
        sleep($sleep);
        $sleep +=30;
        if ($sleep >= 300) { $sleep = 30; }
        }
    }
########## END wait_for_available_host


sub wait_stable_file {
    my ($i, $j, $k);
    my $f = shift;
    return if ($stable_files{$f});
    return unless (-e $f);
    
    # 看不懂
    if (-e "$f.done") 
    { 
        $stable_files{$f} = 1; 
        return; 
    }
    
    my $size0 = -s $f;
    while(1) {
        sleep(10);
        my $size1 = -s $f;
        if ($size0 == $size1) 
        { 
            $stable_files{$f} = 1; 
            last; 
        }
        else 
        {
            $size0 = $size1; 
        }
    }
}
########## END wait_stable_file


# 刷新restart文件
sub write_restart {
    my ($i, $j, $k);
    open(RES, "> $restart_file") || die;
    
    for ($i=0; $i<$command_no; $i++) {
        print RES "$commands[$i]\n$command_status[$i]\n";
    }
    close(RES);
}
########## END write_restart


# 分配执行命令
sub assign_commands {
    my ($i, $j, $k);
    my $cmd;
    my ($idb, $idbo, $jdb, $idbout, $idblog);
    
    $command_no = 0;
    # 切分fasta输入文件
    # cd-hit-div -i 292S_oral_allgene.fa
    #            -o 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div
    #            -div 64

    $cmd = "$cd_hit_div_exe -i $in -o $indiv -div $seg_no";
    
    # 添加切分命令到commands数组里面
    push(@commands, $cmd);
    # 给该命令分配运行状态
    push(@command_status, "wait");
    # 对要运行的命令计数
    $command_no++;

    # 根据切分次数创建子目录，存放中间过程的比对结果
    # [0, 63)
    # mkdir 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div-0-
    
    # mdkir 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div-1-
    # ...
    # mkdir 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div-62-
    # mdkir 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div-63-
    
    for ($i = 0; $i < $seg_no; $i++) 
    {
        $idb    = "$indiv-$i";       
        $idblog = "$indiv-$i-/log";
        `mkdir $indiv-$i-` unless (-e "$indiv-$i-");
        # compare to previous segs
        for ($j = 0; $j < $i; $j++) 
        {
            ## 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div-0-/o
            ## 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div-0-/vs.0
            ## 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div-0-/vs.1
            ## 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div-0-/vs.2
            ## ......
            ## 292S_oral_nr.cd-hit-para-tmp/292S_oral_allgene.fa.div-0-/vs.(i-1)

            $jdb = "$indiv-$j-/o";      # 输入
            $idbo = "$indiv-$i-/vs.$j"; # 输出
            
            # i = 0 nothing
            
            # i = 1 all.div-0/o vs all.div-1 -> all.div-1-/vs.0
            
            # i = 2 all.div-0/o vs all.div-2 -> all.div-2-/vs.0
            #       all.div-1/o vs all.div-2 -> all.div-2-/vs.1
            
            # i = 3 all.div-0/o vs all.div-3 -> all.div-3-/vs.0
            #       all.div-1/o vs all.div-3 -> all.div-3-/vs.1
            #       all.div-2/o vs all.div-3 -> all.div-3-/vs.2

            # i = i all.div-0/o vs all.div-i -> all.div-i-/vs.0
            #       all.div-0/o vs all.div-i -> all.div-i-/vs.1
            #       ......
            #       all.div-0/o vs all.div-i -> all.div-i-/vs.i

            $cmd = "$cd_hit_2d_exe -i $jdb -i2 $idb -o $idbo $arg_pass >> $idblog";
            push(@commands，$cmd);
            push(@command_status, "wait");
            $command_no++;
            $idb = $idbo;
        }
        # self comparing
        $cmd = "$cd_hit_exe -i $idb -o $indiv-$i-/o $arg_pass >> $idblog";
        push(@commands, $cmd);
        push(@command_status, "wait");
        $command_no++;
    }
}

# 当这个cd-hit-para.v2.pl意外退出或者有些正在run的程序出现问题时，可以手动把run改为wait
# 这样相当于手动刷新了restart文件
# 然后再执行ch-hit-para.v2.pl，指定restart文件
sub read_restart {
    $command_no = 0;
    open(RRRR, "$restart_in") || die;
    my $ll;
    while ($ll = <RRRR>) {
        chop($ll);
        push(@commands, $ll);
        $ll = <RRRR>;
        chop($ll);
        push(@command_status, $ll);
        $command_no++;
    }
    close(RRRR);
}
########## END read_restart


sub print_usage {
    print <<EOD;
Usage: $script_name options
        This script divide a big clustering job into pieces and submit
        jobs to remote computers over a network to make it parallel. 
        After all the jobs finished, the script merge the clustering
        results as if you just run a single cd-hit or cd-hit-est.

        You can also use it to divide big jobs on a single computer if
        your computer does not have enough RAM (with -L option).

Requirements:
      1 When run this script over a network, the directory where you 
        run the scripts and the input files must be available on 
        all the remote hosts with identical path.
      2 If you choose "ssh" to submit jobs, you have to have 
        passwordless ssh to any remote host, see ssh manual to
        know how to set up passwordless ssh.
      3 I suggest to use queuing system instead of ssh, 
        I currently support PBS and SGE
      4 cd-hit cd-hit-2d cd-hit-est cd-hit-est-2d
        cd-hit-div cd-hit-div.pl must be in same directory where
        this script is in.
Options 

     -i input filename in fasta format, required
     -o output filename, required
    --P program, "cd-hit" or "cd-hit-est", default "cd-hit"
    --B filename of list of hosts, 
        requred unless -Q or -L option is supplied
    --L number of cpus on local computer, default $local_cpu
        when you are not running it over a cluster, you can use 
        this option to divide a big clustering jobs into small
        pieces, I suggest you just use "--L 1" unless you have
        enough RAM for each cpu
    --S Number of segments to split input DB into, default $seg_no
    --Q number of jobs to submit to queue queuing system, default $queue
        by default, the program use ssh mode to submit remote jobs
    --T type of queuing system, "PBS", "SGE" are supported, default $queue_type
    --R restart file, used after a crash of run
     -q qsub parameter [-q wc_queue_list] default all availabile queues
     -P qsub parameter [-P project_name] default not
    --l qsub parameter [-l resource_list]: vf=xxG[,p=xx,...] (default vf=1G)
     -h print this help

More cd-hit/cd-hit-est options can be speicified in command line

    Questions, bugs, contact Weizhong Li at liwz\@sdsc.edu

EOD
    exit;
}
#### END print_usage
