use warnings;
use strict;
 
my $ori = "00.data";# data file you want to trim or modify 
my $ref = "shearstrain.data";#the ref data for modifying $org

my $retype_natom = `grep "atoms" ./$ref|awk '{print \$1}'`;
chomp $retype_natom;
print "$retype_natom\n";
my @retypeID = `grep -v '^[[:space:]]*\$' ./$ref|grep -A $retype_natom "Atoms"|grep -v "Atoms"|grep -v -- '--'|awk '{print \$1}'`;
map { s/^\s+|\s+$//g; } @retypeID; 

my $natom = `grep "atoms" ./$ori|awk '{print \$1}'`;
chomp $natom;
print "natom\n";
my @ID = `grep -v '^[[:space:]]*\$' ./$ori|grep -A $natom "Atoms"|grep -v "Atoms"|grep -v -- '--'`;
map { s/^\s+|\s+$//g; } @ID;
my $count = 0;

for (0..$#ID){
    my $coors = $ID[$_];
    my @temp = split(/\s+/,$coors);
    if($temp[0] ~~ @retypeID){
        print "$coors\n";
        $temp[1] = $temp[1] + 2;#type
        $count++;
        my $modifed = join (" ",@temp);
        $ID[$_] = "$modifed";
    }
}

for my $m (@ID){

    print "$m\n";
}
#my $here_doc =<<"END_MESSAGE";
## $in
#
#$nat atoms
#$atom_types atom types
#
#$lmp_cell
#
#Masses
#
#$mass4data
#
#Atoms  # atomic
#
#$coords4data
#END_MESSAGE
#
#open(FH, "> $path/$filename.data") or die $!;
#print FH $here_doc;
#close(FH);