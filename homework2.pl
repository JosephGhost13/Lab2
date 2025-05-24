
#!/usr/bin/perl
use strict;
use warnings;

# Πίνακας γενετικού κώδικα
my %codon_table = (
    'ATA'=>'I', 'ATC'=>'I', 'ATT'=>'I', 'ATG'=>'M',
    'ACA'=>'T', 'ACC'=>'T', 'ACG'=>'T', 'ACT'=>'T',
    'AAC'=>'N', 'AAT'=>'N', 'AAA'=>'K', 'AAG'=>'K',
    'AGC'=>'S', 'AGT'=>'S', 'AGA'=>'R', 'AGG'=>'R',
    'CTA'=>'L', 'CTC'=>'L', 'CTG'=>'L', 'CTT'=>'L',
    'CCA'=>'P', 'CCC'=>'P', 'CCG'=>'P', 'CCT'=>'P',
    'CAC'=>'H', 'CAT'=>'H', 'CAA'=>'Q', 'CAG'=>'Q',
    'CGA'=>'R', 'CGC'=>'R', 'CGG'=>'R', 'CGT'=>'R',
    'GTA'=>'V', 'GTC'=>'V', 'GTG'=>'V', 'GTT'=>'V',
    'GCA'=>'A', 'GCC'=>'A', 'GCG'=>'A', 'GCT'=>'A',
    'GAC'=>'D', 'GAT'=>'D', 'GAA'=>'E', 'GAG'=>'E',
    'GGA'=>'G', 'GGC'=>'G', 'GGG'=>'G', 'GGT'=>'G',
    'TCA'=>'S', 'TCC'=>'S', 'TCG'=>'S', 'TCT'=>'S',
    'TTC'=>'F', 'TTT'=>'F', 'TTA'=>'L', 'TTG'=>'L',
    'TAC'=>'Y', 'TAT'=>'Y', 'TAA'=>'_', 'TAG'=>'_', 'TGA'=>'_'
);

# Διαβάζει την αλληλουχία από το χρήστη ή αρχείο
my $sequence = join("", <>);
$sequence =~ s/\s+//g;
$sequence = uc($sequence);

# Υπολογισμός reverse complement
sub reverse_complement {
    my ($seq) = @_;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return reverse $seq;
}

# Εύρεση και μετάφραση ORFs
sub find_orfs {
    my ($dna, $strand) = @_;
    my @orfs;

    for my $frame (0..2) {
        for (my $i = $frame; $i < length($dna) - 2; $i += 3) {
            my $codon = substr($dna, $i, 3);
            if ($codon eq 'ATG') {
                for (my $j = $i + 3; $j < length($dna) - 2; $j += 3) {
                    my $stop = substr($dna, $j, 3);
                    if ($stop eq 'TAA' or $stop eq 'TAG' or $stop eq 'TGA') {
                        my $orf = substr($dna, $i, $j + 3 - $i);
                        my $protein = translate($orf);
                        push @orfs, { start => $i+1, end => $j+3, strand => $strand, protein => $protein };
                        last;
                    }
                }
            }
        }
    }
    return @orfs;
}

# Μετάφραση DNA σε πρωτεΐνη
sub translate {
    my ($seq) = @_;
    my $protein = "";
    for (my $i = 0; $i < length($seq) - 2; $i += 3) {
        my $codon = substr($seq, $i, 3);
        $protein .= $codon_table{$codon} // 'X';
    }
    return $protein;
}

# Εκτέλεση
my @forward_orfs = find_orfs($sequence, 'forward');
my $rev_comp = reverse_complement($sequence);
my @reverse_orfs = find_orfs($rev_comp, 'reverse');

# Εκτύπωση ORFs
foreach my $orf (@forward_orfs, @reverse_orfs) {
    print "Strand: $orf->{strand}, Start: $orf->{start}, End: $orf->{end}\n";
    print "Protein: $orf->{protein}\n\n";
}
