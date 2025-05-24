
use strict;
use warnings;

while (<>) {
    if (/^[ATCG]+$/i) {

        print "Valid DNA sequence\n";

    }
    else {
        print "Invalid DNA sequence\n";
    }


}