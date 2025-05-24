

 my $protein_line = "MNVEHE _123! LLVEE \$"; 
 my $protein_line =~ s/[^A-Z]//g;
print "Cleaned sequence: $protein_line\n";