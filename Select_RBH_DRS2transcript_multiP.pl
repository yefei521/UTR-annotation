use threads; 
use File::Spec; 
use threads::shared; 
use Fcntl qw(:flock); 
my$dir=$ENV{'PWD'}; 
open(IN2,"<","$dir/Manual_check-total-gene.gff3_Right_NonCodingRNA");
open(IN4,"<","$dir/Manual_check-total-gene.gff3_Right_UTR.gff3");
foreach(<IN4>){
   if(/^(chr_\d+)\s+.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?)\s/){#chr_088 GSAman  mRNA    461906  463746  .       -       .       Parent=YF00013595;ID=YF00013595.t1
    $Mangene{$5}="$1\t$2\t$3\t$4"; #print "$5\t$Mangene{$5}\n";
   }elsif(/^(chr_\d+)\s+.+?CDS\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?Parent=(.+?);/){
    push@{$CDS{"Man"}{$5}},$2,$3;
   }elsif(/^(chr_\d+)\s+.+?\s(.+?UTR)\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?Parent=(.+?);/){
    push@{$UTR{$6}{$2}},$_;
   }
}
$n=keys%Mangene;
print "Man total gene:$n\n";

foreach(<IN2>){
   if(/^(chr_\d+)\s+.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?);/){
       $Anno_NonCod{$5}.=$_;
   }elsif(/^(chr_\d+)\s+.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?Parent=(.+?);/){
      $Anno_NonCod{$5}.=$_;
   }
}
$n=keys%Anno_NonCod;
print "Total annotated noncoding RNA:$n\n";

my @file :shared=glob"$dir/temp_200p/*"; 
foreach my $i ( 0 .. 199){     
    $p=$i;   
    print "$p\n";  
    if( fork() ){        wait if($i+1 >= 20);     
    }else{     
       print "START thread: $p ......\n"; 
       &MATCH($p,$dir);        
        #flock OUT, LOCK_EX;        
        #print OUT"$p\t$uniq\t$total\n";        
        #flock OUT, LOCK_UN;        
        exit();     
        }     
    sleep 0.5;
}  
wait;  
while( wait != -1 ){     sleep 5; } 

sub MATCH{    
    my($zmwp,$dir)=@_;  $name=$zmwp;
    open(my$zmwp,"<$dir/temp_200p/$zmwp") or die;   
    while(<$zmwp>){
       next if $_=~/#/;
       @temp=split(/\s+/,$_);
       next if $temp[10]>1e-3;
       if(!exists$DRS2Man{$temp[0]}){$DRS2Man{$temp[0]}{$temp[1]}++;}
       push@{$DRS2Man_DRS_region{$temp[0]}{$temp[1]}},$temp[6],$temp[7];
       push@{$DRS2Man_man_region{$temp[0]}{$temp[1]}},$temp[8],$temp[9];
    }
    open(OUT1,">","$dir/temp_200p/${name}_nonCodingRNA") or die;
    open(OUT2,">","$dir/temp_200p/${name}_AS") or die;
    open(OUT3,">","$dir/temp_200p/${name}_AS.bed") or die;
    open(NCR_Long,">","$dir/temp_200p/${name}_unique_Longest_nonCodingRNA") or die;
    foreach $drs(keys%DRS2Man_DRS_region){
      foreach $man(keys%{$DRS2Man_DRS_region{$drs}}){
        next if !exists$DRS2Man{$drs};
        @DRS_region=@{$DRS2Man_DRS_region{$drs}{$man}};
        @Man_region=@{$DRS2Man_man_region{$drs}{$man}};
    #   print "$drs:@DRS_region\n$man:@Man_region\n";
        if($Man_region[0]>$Man_region[1]){
           if(keys%{$DRS2Man_DRS_region{$drs}}==1){ $NonCoding{$man}{$drs}++; }
           print OUT1"$drs|$man\t$Mangene{$man}\n" or die;
           last;
        }
        next if !exists $Mangene{$man};
        foreach $utr(keys%{$UTR{$man}}){
              foreach $l(@{$UTR{$man}{$utr}}){
                 if($l=~/^(chr_\d+)\s+.+?\s(.+?UTR)\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?Parent=($man);/){
                     $UTR_len{$man}{$2}+=$4-$3+1;
                 }
              }
        }
        @sort_Man_region=sort{$a<=>$b}@Man_region;
        @sort_DRS_region=sort{$a<=>$b}@DRS_region;
        $Mangene{$man}=~/chr_\d+\s+(\d+)\s+(\d+)\s+(.)/; $gene_start=$1; $gene_end=$2; $ori=$3;
       # if($ori eq "+"){ $cut1=$sort_DRS_region[0]+($UTR_len{$man}{"five_prime_UTR"}-$sort_Man_region[0]);
       #                  $cut2=$sort_DRS_region[-1]-($gene_end-$gene_start-$sort_Man_region[-1]);
        #}else{ $cut1=$sort_DRS_region[0]+($gene_end-$gene_start-$sort_Man_region[-1]);
        #       $cut2=$sort_DRS_region[-1]-($UTR_len{$man}{"five_prime_UTR"}-$sort_Man_region[0]);
        #}
        $i=0; while($i<@Man_region){$Man_region[$i]=$gene_start+$Man_region[$i]; $i++;  } 
        @sort_Man_region=sort{$a<=>$b}@Man_region;
        @CDS=@{$CDS{"Man"}{$man}};
        if($drs=~/24fb15e4-7b8c-4242-871f-9a9f2ace20e8/){ print "$man\t$drs\ngene_start $gene_start\tgene_end $gene_end\nCDS:@CDS\nMan_region:@Man_region\nsort_Man_region:@sort_Man_region\n"; }
        if($sort_Man_region[0]<$CDS[1] && $sort_Man_region[-1]>$CDS[-2]){
           print  OUT2"$drs|$man\t$Mangene{$man}\n" or die;
           print  OUT3"$drs\t$sort_DRS_region[0]\t$sort_DRS_region[-1]\t$man\t.\t$ori\n" or die;
        }
      }
    }
    foreach$man(keys%NonCoding){
       foreach $drs(keys%{$NonCoding{$man}}){
          @DRS_region=@{$DRS2Man_DRS_region{$drs}{$man}};
          @sort_DRS_region=sort{$a<=>$b}@DRS_region;
          $len=$sort_DRS_region[-1]-$sort_DRS_region[0]+1;
          $Longest{$man}{$drs}=$len;
      }
      #print NCR_Long "$man\t";
      foreach $sorted_Longest_drs(sort{$Longest{$man}{$b}<=>$Longest{$man}{$a}}keys%{$Longest{$man}}){
          print NCR_Long "$man\t$sorted_Longest_drs\t$Longest{$man}{$sorted_Longest_drs}\n";
      }
      #print NCR_Long "\n";
    }
    close OUT1; close OUT2; close OUT3; close NCR_Long;
    return();
}

print "FINISH!>>>>>";
