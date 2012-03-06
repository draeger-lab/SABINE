#!/usr/bin/perl -w

%seenIDs=();

$validID=0;

while (<>) {
  
  if (m/>(\S+)\s+(\S+)/) {
    $key="$1:$2";
    if (!exists($seenIDs{$key})) {
      $seenIDs{$key}=1;
      print ">$1:$2\n";
      $validID=1;
    } else {
      $validID=0;
    }
  } else {
    if ($validID) {
      chop;
      @sequence=split(//,uc($_));
      print "@sequence\n";
    }
  }
}
