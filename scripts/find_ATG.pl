#!/usr/bin/env perl
use strict;
use warnings;
# QQ group: perlchina
# question: find longest ATG+ sequences

my $seq = "ATGATGASFSAGATGATGATGSFAATGATGATGATGDSFS";

my @atg = $seq =~ /((ATG)+)/g;
my @atg_len = map { length($_) } @atg;
print "@atg\n";
print "@atg_len\n\n";

print((sort {$b cmp $a} ($seq =~ /(?:ATG)+/g))[0]);
print "\n\n";

my @atg_2 = $seq =~ /(?:ATG)+/g;
my @atg_len_2 = map { length($_) } @atg_2;
print "@atg_2\n";
print "@atg_len_2\n";