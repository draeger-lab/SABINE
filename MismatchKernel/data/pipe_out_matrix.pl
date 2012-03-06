#!/usr/bin/perl -w

my $output = "$ARGV[0]";
open(OUT, ">$output");
my $flag = 0;
while(<STDIN>){
   if ($flag == 0){
      $flag = 1;
   }
   else{
    my $line = $_;
    my @tokens  = split("\t",$line); 
    for (my $i = 1; $i < scalar(@tokens); $i++){
       print OUT "$tokens[$i] ";
    }
   }
}
close(OUT);

