#!/usr/bin/perl -w

%seenIDs=();

while (<>) {
  
  if (m/>(\S+)\s+(\S+)/) {
    $key="$1:$2";
    if (!exists($seenIDs{$key})) {
      $seenIDs{$key}=1;
      print "$key\n";
    }
  }
}

#while (($key,$value)=each(%seenIDs)) {
#  print "$key\n";
#}
