#!/usr/bin/perl -w

%counts=();

#$counts{"A"}=0;
#$counts{"C"}=0;
#$counts{"G"}=0;
#$counts{"T"}=0;

$totalCount=0;
for ($i=0; $i<scalar(@ARGV); $i++) {

  open INPUT, "<$ARGV[$i]";

  while (<INPUT>) {

    chop;
    $currentLine=$_;
    
    @currentSymbols=split(/ /,$currentLine);
    #print "@currentSymbols";

    for ($j=0; $j<scalar(@currentSymbols); $j++) {
      if (!($currentSymbols[%j]=~m/^>/)) {
	#print "$currentSymbols[$j]\n";
	$counts{"$currentSymbols[$j]"}++;
	$totalCount++;
      }
    }
  }
  close INPUT;
}

#print "A\t".$counts{"A"}/$totalCount."\n";
#print "C\t".$counts{"C"}/$totalCount."\n";
#print "G\t".$counts{"G"}/$totalCount."\n";
#print "T\t".$counts{"T"}/$totalCount."\n";


while (($key, $value) = each(%counts)) {
  print "$key\t".$counts{$key}/$totalCount."\n";
}
