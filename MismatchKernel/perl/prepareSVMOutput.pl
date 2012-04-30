#!/usr/bin/perl -w

while (<>) {

  if (m/^\#/) {
  }
  else {
    m/(\S+)\s+\S+\s+(\S+)/;
    print "$1\t$2\n";
  }
}
