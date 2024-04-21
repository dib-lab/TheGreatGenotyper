#!/usr/bin/perl

use strict;
use warnings;

# Process input line by line
while (my $line = <STDIN>) {
    chomp $line;
    
    # If the line is a header line, print it as is
    if ($line =~ /^#/) {
        print "$line\n";
    } else {
        # Split the line into columns on tab character
        my @fields = split(/\t/, $line);
        
        # Replace ';' with '_' in the ID field (third column)
        $fields[2] =~ s/;/_/g;

        # Join the fields back together and print
        print join("\t", @fields) . "\n";
    }
}
