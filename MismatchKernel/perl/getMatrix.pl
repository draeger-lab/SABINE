#!/usr/bin/perl -w

$beforeMatrix=1;

while (<>) {
  if ($beforeMatrix) {
    if (m/Kernel Matrix/) {
      $beforeMatrix=0;
    }
  } else {
    print;
  }
}
