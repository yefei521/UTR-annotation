open(IN1,"<","DRS2Manual_otf7") or die;
my$dir=$ENV{'PWD'};

$cmd=`mkdir temp_200p`;
print $cmd;
foreach my$i(0 .. 200){
   open($i,">","$dir/temp_200p/$i") or die;
}

$line=22417986;
while(<IN1>){
    next if $_=~/#/;
    $n++; $p=int(($n*200)/$line);
    print $p "$_";  
}
foreach my$i(0 .. $p){
   close $i;
}

