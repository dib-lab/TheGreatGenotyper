#!/usr/bin/perl
use strict;
use warnings;

# Flag to indicate when we are in the header section
my $in_header = 1;

while (my $line = <STDIN>) {
    # Check if we're still in the header
    if ($in_header) {
        # Output the line as is
        print $line;

        # Add new header lines after the file format line
        if ($line =~ /^##fileformat=/) {
            print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n";
        }

        # Check if we have reached the end of the header
        if ($line =~ /^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO/) {
            $in_header = 0;
        }
    } else {
        # If we're past the header, output the line as is
        print $line;
    }
}
