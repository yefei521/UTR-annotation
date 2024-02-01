$dir=$ENV{'PWD'};
open(IN1,"<","$dir/WT_mix_seqkit2DNA_MaxIntron2k.sam");
#open(IN1,"<","$dir/temp.sam");
#open(IN2,"<","$dir/WT_mix_seqkit2DNA_MaxIntron2k_sort.bed");
open(IN2,"<","$dir/Manual_check-total-gene.gff3_Right_UTR_add-loqQ-1271.gff3");

@IN2=<IN2>; close IN2;
foreach(@IN2){
    if(/^(chr_\d+)\s+.+?\s+gene\s+(\d+)\s+(\d+)\s.+?\s+(.)\s.+?ID=(.+?);/){
       #push@{$gene{$5}},$1,$4,$2,$3;
       push@{$gene{$1}{$5}},$4,$2,$3; $n++;
    }
}
print "$n\n";

while(<IN1>){
   @temp=split(/\t/,$_);
   if($temp[1]==0){
       #print "$temp[0]\n";
       $match=$temp[5]; $pos=$temp[3];  $strand{$temp[0]}="+"; 
       while($match=~/^(\d+)([MDSIHN])(.+)/){  #6167d4cf-9149-4179-8de5-258a6d61c8b5    0       chr_133 432544  60      27M2I4M1I29M5D5M1D2M1D24M1I48M2D53M2I43M1D21M20S        *       0       0       TGCCCTTATGCCC
          $match=$3; $MM=$2; $l=$1;  
          if($MM=~/N|M|D/){ $intron_5=$pos; $pos+=$l; if($MM=~/N/){ $intron_3=$pos; push@{$DRS{$temp[0]}{"Intron"}{$temp[2]}},"$intron_5-$intron_3"; } }  
       }
       #print "$temp[0]\t$pos\n";
       $mid=int(($temp[3]+$pos)/2); $mark=0;
       foreach$name(keys%{$gene{$temp[2]}}){ 
            @gene=@{$gene{$temp[2]}{$name}};
            if($gene[0] eq "+" && $mid>= $gene[1] && $mid<=$gene[2]){ push@{$gene2reads{$name}},$temp[0]; $mark++;last;  }
            if($gene[0] eq "-" && $mid>= $gene[1] && $mid<=$gene[2]){ push@{$gene2NAT{$name}},$temp[0]; $mark++; last; }
        }
        push@intergenic_read,$temp[0] if $mark==0;
   }elsif($temp[1]==16){
       $match=$temp[5]; $pos=$temp[3];  $strand{$temp[0]}="-";  
       while($match=~/^(\d+)([MDSIHN])(.+)/){  #6167d4cf-9149-4179-8de5-258a6d61c8b5    0       chr_133 432544  60      27M2I4M1I29M5D5M1D2M1D24M1I48M2D53M2I43M1D21M20S        *       0       0       TGCCCTTATGCCC
          $match=$3; $MM=$2; $l=$1;  
          if($MM=~/N|M|D/){ $intron_5=$pos; $pos+=$l; if($MM=~/N/){ $intron_3=$pos; push@{$DRS{$temp[0]}{"Intron"}{$temp[2]}},"$intron_5-$intron_3"; } }               
       }
       $mid=int(($temp[3]+$pos)/2); $mark=0;
       foreach$name(keys%{$gene{$temp[2]}}){ 
            @gene=@{$gene{$temp[2]}{$name}};
            if($gene[0] eq "-" && $mid>= $gene[1] && $mid<=$gene[2]){ push@{$gene2reads{$name}},$temp[0]; $mark++;last;  }
            if($gene[0] eq "+" && $mid>= $gene[1] && $mid<=$gene[2]){ push@{$gene2NAT{$name}},$temp[0]; $mark++; last; }
        }   
        push@intergenic_read,$temp[0] if $mark==0;
   }
}

foreach(@IN2){
   if(/^(chr_\d+)\s+.+?\s+gene\s+(\d+)\s+(\d+)\s.+?\s+(.)\s.+?ID=(.+?);/){ #chr_124 GSAman  gene    190461  192058  .       +       .       ID=YF00019807;Name=YF00019807
       $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5;
       #sense reads counts
       @reads=@{$gene2reads{$name}};
       foreach$id(@reads){ $splice_site{"sense"}{$name}{$splice_site}=0;
           if($strand{$id} eq $ori){
               foreach $splice_site(@{$DRS{$id}{"Intron"}{$chr}}){ $splice_site{"sense"}{$name}{$splice_site}++;  }  
               $countreads{"sense"}{$name}++;
           }
       }
       #antisense reads count 
       @reads=@{$gene2NAT{$name}};
       foreach$id(@reads){ $splice_site{"antisense"}{$name}{$splice_site}=0;
           if($strand{$id} ne $ori){
               foreach $splice_site(@{$DRS{$id}{"Intron"}{$chr}}){ $splice_site{"antisense"}{$name}{$splice_site}++;  }
               $countreads{"antisense"}{$name}++;
           }
       } 
   }
}

open(OUT1,">","$dir/Gene-NAT_intron-reads_count.txt");
open(OUT2,">","$dir/Gene-NAT_transcript-reads_count.txt");
foreach $name(keys%{$splice_site{"sense"}}){
   print OUT1"$name\tsense\t"; $intron_num=0;
   foreach $ss(keys%{$splice_site{"sense"}{$name}}){
      print OUT1"$ss:$splice_site{sense}{$name}{$ss}\t";
      $intron_num++;
   }
   print OUT1"\n";
   print OUT2"$name\tsense\tintron-num:$intron_num\treads-num:$countreads{sense}{$name}\n";
}

foreach $name(keys%{$splice_site{"antisense"}}){
   print OUT1"$name\tantisense\t"; $intron_num=0;
   foreach $ss(keys%{$splice_site{"antisense"}{$name}}){
      print OUT1"$ss:$splice_site{antisense}{$name}{$ss}\t";
      $intron_num++;
   }  
   print OUT1"\n";
   print OUT2"$name\tantisense\tintron-num:$intron_num\treads-num:$countreads{antisense}{$name}\n";
}

foreach $id(@intergenic_read){
   print OUT1"NULL\tintergenic\t$id\n";
}

