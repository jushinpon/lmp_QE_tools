#!/usr/bin/perl
=b
https://docs.lammps.org/Howto_triclinic.html
The parallelepiped has its “origin” at (xlo,ylo,zlo) and 
is defined by 3 edge vectors starting from the origin given by
a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); c = (xz,yz,zhi-zlo). 

cfg format:
ITEM: BOX BOUNDS xy xz yz
xlo_bound xhi_bound xy
ylo_bound yhi_bound xz
zlo_bound zhi_bound yz

data format:
0.000000000000      11.327003370000  xlo xhi
0.000000000000      10.970021530071  ylo yhi
0.000000000000      10.529870936067  zlo zhi
0.507912903230       0.714988108603      -0.046041891158  xy xz yz
=cut
use strict;
use warnings;
use Cwd;
use Data::Dumper;
use JSON::PP;
use Data::Dumper;
use List::Util qw(min max);
use Cwd;
use POSIX;
use Parallel::ForkManager;
use lib '.';#assign pm dir
use elements;#all setting package

###parameters to set first
my $currentPath = getcwd();
my @myelement =  ("B","N");#corresponding to lmp type ids
my @datafile = sort `find $currentPath/initial -name "*.data"`;#find all data files
chomp @datafile;
my @temperature = ("300","500","1000");#temperatures for QE_MD, 0 for scf
my @pressure = ("0");#pressure for vc-md, if 0 for scf
my $calculation = "vc-md";#set temperature and pressure to be 0 0 for scf
my $stepsize = 50;#20 ~ 0.98 fs
my $nstep = 200;#how many steps for md for vc-relax
my $pseudo_dir = "/opt/QEpot/SSSP_efficiency_pseudos/";
####end of setting parameters
###get pot setting here!
my $json;
{
    local $/ = undef;
    open my $fh, '<', '/opt/QEpot/SSSP_efficiency.json' or die "no QE pot file\n";#or precision
    $json = <$fh>;
    close $fh;
}
my $decoded = decode_json($json);
#print "$decoded->{\"B\"}->{filename}\n";
#die;
######   cutoff   #######
my @rho_cutoff;
my @cutoff;
my @pot;
my %used_element;
for (@myelement){#unique
    chomp;
    #print "element: $_\n";
    push @rho_cutoff,$decoded->{$_}->{rho_cutoff};
    push @cutoff,$decoded->{$_}->{cutoff};
 #density (g/cm3), arrangement, mass, lat a , lat c
    #my @temp = 
    @{$used_element{$_}} = &elements::eleObj("$_");
    die "no information of element $_ in elements.pm\n" unless (@{$used_element{$_}});
    die "no pot file of element $_\n" unless ($decoded->{"$_"}->{filename});
    push @pot,[$_,${$used_element{$_}}[2],$decoded->{"$_"}->{filename}];
}
# for keeping the largest ones only
my @rho_cutoff1 = sort {$b<=>$a} @rho_cutoff;
my @cutoff1 = sort {$b<=>$a} @cutoff;
my $ntyp = @pot; #element type number
## for here doc
my $rho_cutoff = $rho_cutoff1[0];
my $cutoff = $cutoff1[0];
my $species = "";
my $starting_magnetization = "";
my $counter = 0;
for (@pot){
    $counter++;
    my $temp = join(" ",@{$_});
    #print $temp."\n";
    $species .= "$temp\n";
    $starting_magnetization .= "starting_magnetization($counter) = 0.01\n";
}
chomp ($species,$starting_magnetization);
my $myelement = join ('',@myelement);
`rm -rf disturb`;
`mkdir disturb`;
`mkdir -p disturb/QEinput`;
`rm -rf original`;
`mkdir original`;
for my $id (@datafile){
    my $data_path = `dirname $id`;
    my $data_name = `basename $id`;
    $data_name =~ s/\.data//g;
    chomp ($data_path, $data_name);
    `cp $id ./original/$data_name.lmp`;
    for my $temp (@temperature){
    for my $press (@pressure){
        #print "$id, T,P: $temp,$press\n";
        #`mkdir -p disturb/$data_name`;
         my $output = $data_name ."_T$temp"."_P$press.lmp";
        `atomsk original/$data_name.lmp -disturb 0.05 0.05 0.05 disturb/$output`;
        open my $database ,"< disturb/$output";
      
        my @data =<$database>;
        close $database;
        my %para = (
            xy => "0.0",
            xz => "0.0",
            yz => "0.0",
            natom => "",
            xlo => "",
            xhi => "",
            ylo => "",
            yhi => "",
            zlo => "",
            zhi => "",
            coords => []
        );
        for (@data){
            chomp;
    ####atoms###
            if(/(\d+)\s+atoms/){ 
                $para{natom} = $1;
            }
    ####CELL_PARAMETERS###
            elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+xlo\s+xhi/){
                $para{xlo} = $1;
                $para{xhi} = $2;
            }
            elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+ylo\s+yhi/){
                $para{ylo} = $1;
                $para{yhi} = $2;
            }
            elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+zlo\s+zhi/){
                $para{zlo} = $1;
                $para{zhi} = $2;
            }
            elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+xy\s+xz\s+yz/){
                $para{xy} = $1;
                $para{xz} = $2;
                $para{yz} = $3;
            }
#1 1 4.458517505863 1.201338326940 0.873835074284
            elsif(/\d+\s+(\d+)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)$/){
                my $ele = $myelement[$1-1];
                my $x = $2 - $para{xlo};
                my $y = $3 - $para{ylo};
                my $z = $4 - $para{zlo};
                my $temp = join(" ",($ele, $x, $y, $z));
                #print "$temp\n";
                push @{$para{coords}},$temp;
            }
        }#one data file
        #check all data
        for my $k (sort keys %para){
           die "\$para{$k} is empty for $id\n" unless($para{$k});
        }
        
        my $coords = join("\n",@{$para{coords}});
        #cell parameter
        my $lx = sprintf("%.6f",$para{xhi}-$para{xlo});
        my $ly = sprintf("%.6f",$para{yhi}-$para{ylo});
        my $lz = sprintf("%.6f",$para{zhi}-$para{zlo});
        my $a = join(" ",($lx,"0.00","0.00"));
        my $b = join(" ",(sprintf("%.6f",$para{xy}),$ly,"0.00"));
        my $c = join(" ",(sprintf("%.6f",$para{xz}),sprintf("%.6f",$para{yz}),$lz));
        my $cell = join("\n",($a,$b,$c));
        my $QEoutput = $data_name ."_T$temp"."_P$press.in";
    my %QE_para = (
            calculation => "$calculation",
            output_file => "$currentPath/disturb/QEinput/$QEoutput",
            pseudo_dir => $pseudo_dir,
            coord => $coords,
            temp => $temp,
            press => $press,
            ntyp => $ntyp,
            dt => $stepsize, ###timestep size,
            nat => $para{natom},
            nstep => $nstep,
            cell_para => $cell,
            pot_spec => $species,
            starting_magnetization => $starting_magnetization,
            rho_cutoff => $rho_cutoff,
            cutoff => $cutoff
            );
      
        &QEinput(\%QE_para);
    }#pressure
    }#temperature
}#all data files
   
######here doc for QE input##########
sub QEinput
{

my ($QE_hr) = @_;

my $QEinput = <<"END_MESSAGE";
&CONTROL
calculation = "$QE_hr->{calculation}"
nstep = $QE_hr->{nstep}
etot_conv_thr = 1.0d-5
forc_conv_thr = 1.0d-4
disk_io = '/dev/null'
pseudo_dir = '$QE_hr->{pseudo_dir}'
tprnfor = .true.
tstress = .true.
verbosity = 'high'
dt = $QE_hr->{dt}
/
&SYSTEM
ntyp =  $QE_hr->{ntyp}
occupations = 'smearing'
smearing = 'gaussian'
degauss =   7.3498618000d-03
ecutrho =   $rho_cutoff
ecutwfc =   $cutoff 
ibrav = 0
nat = $QE_hr->{nat}
nosym = .TRUE.
$QE_hr->{starting_magnetization}
nspin = 2
/
&ELECTRONS
conv_thr =   2.0000000000d-10
electron_maxstep = 200
mixing_beta =   4.0000000000d-01
/
&IONS
ion_dynamics = "beeman"
ion_temperature = "rescaling"
tempw = $QE_hr->{temp}
/
&CELL
!press_conv_thr = 0.1
cell_dynamics = "pr"
press = $QE_hr->{press}
cell_dofree = "volume"
/
K_POINTS {automatic}
2 2 2 0 0 0
ATOMIC_SPECIES
$QE_hr->{pot_spec}
ATOMIC_POSITIONS {angstrom}
$QE_hr->{coord}
CELL_PARAMETERS {angstrom}
$QE_hr->{cell_para}
END_MESSAGE

my $file = $QE_hr->{output_file};
open(FH, '>', $QE_hr->{output_file}) or die $!;
print FH $QEinput;
close(FH);
}