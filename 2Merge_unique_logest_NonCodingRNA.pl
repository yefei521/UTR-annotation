$dir=$ENV{'PWD'};
my @file :shared=glob"$dir/temp_200p/*_unique_Longest_nonCodingRNA";
foreach$file(@file){
     open(IN,"<","$file");
     while(<IN>){
       if(/^(.+?)\s+(.+?)\s+(.+?)\s/){
          $Longest{$1}{$2}=$3; 
        }else{ print "$_"; }
     }
}
open(NCR,">","$dir/DRS2Manual_otf7_unique_Longest_nonCodingRNA");
print NCR"gene\tlongestNoncodingRNA_id\tsecond_NoncodingRNA_id\t.....\n";
foreach $man(keys%Longest){
     print NCR"$man\t";
     foreach $sorted_Longest_drs(sort{$Longest{$man}{$b}<=>$Longest{$man}{$a}}keys%{$Longest{$man}}){
          print  NCR"$sorted_Longest_drs\t";
      }
     print NCR "\n";
}

