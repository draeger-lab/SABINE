#!/usr/bin/perl -w

$scalingFactor=(3/2)/log(2);

@columns = ();
%matrix=();

$line=<>;

#$line=~s/\s*//;




@columns=split(/\s+/,$line);

shift(@columns);


@newColumns=("A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", 
"Q", "R", "S", "T", "V", "W", "X", "Y", "Z");

#@backProbs=(.078,.024,.052,.058,.043,.083,.024,.062,.055,.091,.024,.042,.044,.034,.050,.060,.055,.073,.014,.034);
@backProbs=(.0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435, .0435);




#print("Columns @columns\nNew Columns @newColumns\n");


#for ($i=0; $i<scalar(@columns); $i++) {
#  print("$i=$columns[$i]\n");
#}

$currentRow=0;
while (<>) {
  
  @values=split(/\s+/,$_);
  
  shift(@values);

  for ($i=0; $i<scalar(@values); $i++) {
    $key="$columns[$currentRow]"."-"."$columns[$i]";
    $revKey="$columns[$i]"."-"."$columns[$currentRow]";
    #print "$key = $values[$i]\n";
    $matrix{$key}=$values[$i];
    $matrix{$revKey}=$values[$i];
    #print "$matrix{$key}\n";
  }
  #print "\n";
  if (scalar(@values)>0) {
    $currentRow++;
  }
}

for ($i=0; $i<scalar(@newColumns); $i++) {
  $totalProb=0;
  for ($j=0; $j<scalar(@newColumns); $j++) {
    $key="$newColumns[$i]"."-"."$newColumns[$j]";
    $value=$backProbs[$j]*exp($matrix{$key}/$scalingFactor);
    $totalProb+=$value;
    #print "$value ";
    $newMatrix{$key}=$value;
  }
  #print "\nTotal Prob = $totalProb\n";
  for ($j=0; $j<scalar(@newColumns); $j++) {
    $key="$newColumns[$i]"."-"."$newColumns[$j]";
    $newMatrix{$key}=$newMatrix{$key}/$totalProb;
  }
}

#print "\n\n\n";

for ($i=0; $i<scalar(@newColumns); $i++) {
  $totalProb=0;
  $middlekey="$newColumns[$i]"."-"."$newColumns[$i]";
  
  for ($j=0; $j<scalar(@newColumns); $j++) {
    #if ($i==$j) {
    #  print "0 ";
    #} else {
      $key="$newColumns[$i]"."-"."$newColumns[$j]";
      #$value=$newMatrix{$key}/(1-$newMatrix{$middlekey});
      $value=$newMatrix{$key};
      printf("%1.3f ",$value);
      $totalProb+=$value;
    #}
  }

  print "\n";
  #print "Total Prob = $totalProb\n";
}
