use strict;
use warnings;

my $sequence = "ATGCGATGGTATG";
 my $count = () = $sequence =~ /ATG/g;
print "Start codons found : $count\n";
