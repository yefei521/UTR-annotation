$dir=$ENV{'PWD'};
open(IN1,"<","$dir/gffcmp.annotated.gtf");
open(IN4,"<","$dir/Manual_check-total-gene.gff3_Right_UTR.gff3");
foreach(<IN4>){
   if(/^(chr_\d+)\s+.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?)\s/){#chr_088 GSAman  mRNA    461906  463746  .       -       .       Parent=YF00013595;ID=YF00013595.t1
    $Mangene{$5}="$1\t$2\t$3\t$4"; #print "$5\t$Mangene{$5}\n";
    $Mangene{$5}{"gene_start"}=$2; $Mangene{$5}{"gene_end"}=$3;
   }elsif(/^(chr_\d+)\s+.+?exon\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?Parent=(.+?);/){
    push@{$CDS{"Man"}{$5}},$2,$3;
   }elsif(/^(chr_\d+)\s+.+?\s(.+?UTR)\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?Parent=(.+?);/){
    push@{$UTR{$6}{$2}},$_;
   }
}
$n=keys%Mangene;
print "Man total gene:$n\n";

foreach $l(<IN1>){
   if($l=~/^(chr_\d+)\s+.+?transcript\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?transcript_id "(.+?)"; gene_id.+?gene_name "(.+?)";.+?cmp_ref "\6.t1"; class_code "(.)";/){
        if($7 ne "="){ $ref=$6.".t1";  $name=$5; }
   }elsif($l=~/^(chr_\d+)\s+.+?exon\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?transcript_id "($name)"; gene_id/){ #chr_001 StringTie       exon    2024    2780    .       +       .       transcript_id "AS.2.1"; gene_id "AS.2"; exon_number "1";
          push@{$Man2AS{$ref}{$name}},$2,$3;
   }
}
open(OUT1,">","$dir/DRS2Manual_otf7_AS2TGDnome_sorted.bam_AS.gff3");
foreach$man(keys%Man2AS){
   $n=1;
   foreach $as(keys%{$Man2AS{$man}}){
      @manCDS=@{$CDS{"Man"}{$man}}; #shift@manCDS; pop@manCDS;
      @asCDS=@{$Man2AS{$man}{$as}}; #shift@asCDS; pop@asCDS;
      if($manCDS[1]<$asCDS[1]){ $asCDS[0]=$Mangene{$man}{gene_start};    }
      if($manCDS[-2]>$asCDS[-2]){ $asCDS[-1]=$Mangene{$man}{gene_end};    }
      next if @asCDS<=2;
      $n++; 
      if($Mangene{$man}=~/^(chr_\d+)\s+(\d+)\s+(\d+)\s+(\+|-)/){
          $chr=$1; $gene_start=$2; $gene_end=$3; $ori=$4;
          $man=~/^(.+?).t1/; $name=$1;
          #print OUT1 "$chr\tGSAman\tgene\t$gene_start\t$gene_end\t.\t$ori\t.\tID=$man.$n;Name=$man;ManInfo=AS\n";
           print OUT1 "$chr\tGSAman\tmRNA\t$gene_start\t$gene_end\t.\t$ori\t.\tParent=$name.t$n;ID=$name.t$n\n";
           $i=0; while ($i<@asCDS) {   print OUT1 "$chr\tGSAman\texon\t$asCDS[$i]\t$asCDS[$i+1]\t.\t$ori\t.\tParent=$name.t$n;ID=$name.t$n.exon\n"; $i+=2; }
      }else{ print "$man\n";}
   }   
}



