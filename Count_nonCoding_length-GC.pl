my$dir=$ENV{'PWD'};
open(IN1,"<","$dir/");
open(IN2,"<","$dir/");

foreach(<IN1>){
   if(/^(.+?)\s+(.+?)\s+(.+?)\s/){
      next if $2<30;
      $name=$1; $len=$2; $GC=$3;
      $len{$2}++; $GC{$GC}++;
    }
}

foreach$len(sort{$a<=>$b}keys$len){
   if($len<500){ $a{"0-500"}++;
   }elsif($len<1000){ $a{"500-1000"}++;
   }elsif($len<2000){ $a{"1000-2000"}++;
   }elsif($len<3000){ $a{"2000-3000"}++;
   }elsif($len<4000){ $a{"3000-4000"}++;
   }elsif($len<5000){ $a{"4000-5000"}++;
   }elsif($len<6000){ $a{"5000-6000"}++;
   }elsif($len<7000){ $a{"6000-7000"}++;
   }elsif($len<8000){ $a{"7000-8000"}++;
   }elsif($len<9000){ $a{"7000-9000"}++;
   }elsif($len<10000){ $a{"9000-10000"}++;
   }else{ $a{"10000+"}++;
   } 
}

foreach$GC(sort{$a<=>$b}keys$GC){
    
}
