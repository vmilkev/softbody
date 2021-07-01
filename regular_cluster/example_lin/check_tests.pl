#!/usr/bin/env perl
use strict;
use warnings;
   
# FileSearch.pl: Search for lines containing a search word.
(my $filename, my $word) = @ARGV;   # Get filename & search word.

# Create a filehandle called FILE and connect to the file.
open(FILE, $filename) or die "Can't open $filename: $!";

my $notsuccess = "There is an ERROR in ".$filename;
my $success = 0;
  
while (<FILE>) {           # Set $_ to each line of the file
	# print if /\b$word\b/i;  # Match $_ with word, case insensitive

    #if ($_ == $word){
    if (/\b$word\b/i){
        print "$word\n";
        $success = 1;
    }
}

if ($success == 0){
    print "$notsuccess\n";
} 
