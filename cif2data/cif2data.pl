use strict;
use warnings;
use Cwd;
use POSIX;
#use lib '.';
#use elements;

my $currentPath = getcwd();
my @ciffiles = `find $currentPath/binary/BN -name "*.cif"`;
chomp @ciffiles;

`rm -rf data`;
`mkdir data`;

for my $cif (@ciffiles){
    my $data_path = `dirname $cif`;
    my $data_name = `basename $cif`;
    $data_name =~ s/\.cif//g;
    chomp ($data_path, $data_name);

    my $output = "data/$data_name.lmp";
    my $outputdata = "data/$data_name.data";
   # print "$cif: $output\n";
    unlink "./test.lmp";
    system("atomsk $cif -alignx -unskew $output");
    system("mv $output $outputdata")
}