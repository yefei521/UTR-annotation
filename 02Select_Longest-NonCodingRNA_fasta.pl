$dir=$ENV{'PWD'};
open(IN2,"<","$dir/WT_mix.fa");

print "START______\n";
print "Multi threads_____\n";
foreach my $i ( 0 .. 199){
    $p=$i;
    if( fork() ){        wait if($i+1 >= 60);
    }else{
       print "START thread: $p ......\n";
       &MATCH($p,$dir);
       exit();
    }
    sleep 0.5;
}
wait;
while( wait != -1 ){     sleep 5; }

sub MATCH{
    my($zmwp,$dir)=@_;
    print "ppp:$zmwp\n"; 
    $out=$zmwp."out";
    open($zmwp,"<","$dir/temp_200p/${zmwp}_unique_Longest_nonCodingRNA") or die;
    open($out,">$dir/temp_200p/${zmwp}_unique_Longest_nonCodingRNA.fasta") or die;
    while(<$zmwp>){
       if($_=~/^(.+?)\s+(.+?)\s/){
           $name=$1; 
           $fasta{$name}++;
       }else{print "$_";}
    }
    $n=keys%fasta; print "noncRNA number:$n\n";
    $in=$zmwp."in";
    open($in,"<","$dir/WT_mix.fa");
    while(<$in>){
      if(/>(.+?)\s/ && exists$fasta{$1}){
         $name=$1; $mark=1; print $out ">$name\n";
      }elsif(/>(.+?)\s/ && !exists$fasta{$1}){
        $mark=0; 
      }else{  
        if($mark==1){print $out "$_"; }
      }
    }
    return;
}


