
while (<>) {
    if (/^>\S+\|(\w+)\|/) {
    print "Found accession: $1\n";
    }
}