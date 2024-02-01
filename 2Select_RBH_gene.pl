$dir=$ENV{'PWD'};
open(IN1,"<","$dir/updt1ToManual_otf7");
open(IN2,"<","$dir/ManualToupdt1_otf7");
open(IN3,"<","$dir/2-upd-Genome-GFF3-latest-2.gff3_rmtRNA.gtf");
open(IN4,"<","$dir/Manual_check-total-gene.gff3_Right_UTR.gff3");
@IN3=<IN3>; @IN4=<IN4>;
foreach(<IN1>){
   next if $_=~/#/;
   @temp=split(/\s+/,$_);
   next if $temp[10]>1e-3;
   $upd2Man{$temp[0]}++;
   $upd2Man_score{$temp[0]}{$temp[1]}+=$temp[11];
   $upd2Man_length{$temp[0]}{$temp[1]}+=$temp[3]-$temp[4];
}
$n=keys%upd2Man_score;
print "$temp[0]\n";
print "upd2man:$n\n";

foreach(<IN2>){
   @temp=split(/\s+/,$_);
   next if $_=~/#/;
   next if $temp[10]>1e-3;
   $Man2upd{$temp[0]}++;
   $Man2upd_score{$temp[0]}{$temp[1]}+=$temp[11];
   $Man2upd_length{$temp[0]}{$temp[1]}+=$temp[3]-$temp[4];
}
$n=keys%Man2upd_score;
print "man2upd:$n\n";

foreach(@IN3){
   #print "$_";
   if(/^(chr_\d+)\s+.+?transcript\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?transcript_id "(.+?.t1)";/){ #chr_001 AUGUSTUS        transcript      226     854     .       +       .       transcript_id "g18410.t1"; gene_id "g18410";
      $updgene{$5}="$1\t$2\t$3\t$4"; #print "$5\t$updgene{$5}\n";
   }elsif(/^(chr_\d+)\s+.+?CDS\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?transcript_id "(.+?.t1)";/){
    push@{$CDS{"upd"}{$5}},$2,$3;
   }
}

$n=keys%updgene;
print "$5\n";
print "upd total gene:$n\n";

foreach(@IN4){
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

foreach $upd(keys%updgene){
    next if exists$upd2Man{$upd};
    push@{$type{"Delete"}{$upd}},$upd;
    $delete++;
}

foreach $man(keys%Mangene){
    next if exists$Man2upd{$man};
    push@{$type{"Denovo"}{$man}},$man;
    $denovo++;
}


foreach$upd(keys%upd2Man_score){
     $updgene{$upd}=~/^(.+)\s+(.+)\s+(.+)\s+(.+)$/; $upd_chr=$1; $upd_start=$2; $upd_end=$3; $upd_ori=$4;
     undef@temp;undef@part;
     foreach$man(keys%{$upd2Man_score{$upd}}){
         $Mangene{$man}=~/^(.+)\s+(.+)\s+(.+)\s+(.+)$/; $man_chr=$1; $man_start=$2; $man_end=$3; $man_ori=$4;
         #print "$upd\t$man\t$upd2Man_score{$upd}{$man}\n";
        if($upd_chr eq $man_chr ){
              $mark=0;
              if($man_start>$upd_start && $man_start<$upd_end ){$mark++;}
              if($man_end>$upd_start && $man_end<$upd_end ){$mark++;}
              if($mark>0){ push@{$temp{$upd}},$man; push@part,$man }
        }
     }
     if(@{$temp{$upd}}>1){ @{$type{"Split"}{$upd}}=@{$temp{$upd}}; $split++; push@forExon_Alt,$upd; push@splited,@part; }
}
$splited_parent=@forExon_Alt;
$splited_child=@splited;
undef@part;

foreach$man(keys%Man2upd_score){
     $Mangene{$man}=~/^(.+)\s+(.+)\s+(.+)\s+(.+)$/; $man_chr=$1; $man_start=$2; $man_end=$3; $man_ori=$4;
     undef@temp; undef@part;
     foreach$upd(keys%{$Man2upd_score{$man}}){
         $updgene{$upd}=~/^(.+)\s+(.+)\s+(.+)\s+(.+)$/; $upd_chr=$1; $upd_start=$2; $upd_end=$3; $upd_ori=$4;
         #print "$man\t$upd\t$Man2upd_score{$man}{$upd}\n";
        if($upd_chr eq $man_chr ){
              $mark=0;
              if($upd_start>$man_start && $upd_start<$man_end ){$mark++;}
              if($upd_end>$man_start && $upd_end<$man_end ){$mark++;}
              if($mark>0){ push@{$temp{$man}},$upd; push@part,$upd; }
        }
     }
     if(@{$temp{$man}}>1){ @{$type{"Merge"}{$man}}=@{$temp{$upd}}; $merge++;  push@forExon_Alt,@part;  }
}

$n=@forExon_Alt-$splited_parent;
print "delete:$delete\ndenovo:$denovo\n";
print "split:$split split2 $splited_child\nmerge:$n merged2 $merge\n";

open(OUT1,">","$dir/Exon_alternation_gene.txt");
foreach$upd(keys%upd2Man){
     next if grep(/$upd/,@forExon_Alt);
     undef@man;
     foreach$man(keys%{$upd2Man_score{$upd}}){ 
         push@man,$man;
     }
     next if @man>1;
     $upd2man_one++;
     $man=$man[0];
     @cds_upd=@{$CDS{"upd"}{$upd}}; 
     if(@cds_upd>2){pop@cds_upd; shift@cds_upd;}
     @cds_man=@{$CDS{"Man"}{$man}}; 
     if(@cds_man>2){pop@cds_man; shift@cds_man;}
     if(@cds_upd~~@cds_man){  
     }else{ $Exon_Alt++; print OUT1 "$upd|$man\n@cds_upd\n$man|$upd\n@cds_man\n"; }
}

print "one2one num:$Exon_Alt\n";
print "Exon alternation num:$Exon_Alt\n";

##----------------UTR num
foreach $man(keys%UTR){
    $mark5=0; $mark3=0;
    foreach $line(keys%{$UTR{$man}}){
      if(grep(/five/,@{$UTR{$man}{$line}})){  $mark5++;}
      if(grep(/three/,@{$UTR{$man}{$line}})){ $mark3++;}
    }
    if($mark5>0 && $mark3>0){$two_utr++; 
    }elsif($mark5>0 && $mark3==0){ $utr5++;
    }elsif($mark5==0 && $mark3>0){ $utr3++; 
    }elsif($mark5==0 && $mark3==0){ $none++; }
}

print "gene with UTR5+UTR3 num:$two_utr\n";
print "gene only with UTR5 num:$utr5\n";
print "gene only with UTR3 num:$utr3\n";
print "gene has not UTR5+UTR3 num:$none\n";











