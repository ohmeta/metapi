#/usr/bin/env perl
use strict;
use warnings;
my $a = "hello/world";
my $b = join('', $a, "/summary");
print "$a\n";
print "$b\n";
print "$ARGV[0]\n";
print "$ARGV[1]\n";
print "$ARGV[2]\n";
print "$ARGV[3]\n";

# eg : perl perl_test.pl a b c d
# output:
# hello/world
# hello/world/summary
# a
# b
# c
# d