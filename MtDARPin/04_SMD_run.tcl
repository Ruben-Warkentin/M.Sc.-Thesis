#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# MtDARPin vs  96 T3Q molecules, arranged around it. -- SMD

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################
set xscFile [glob "MtDARPin\_\_T3Q\_*\_\_solvated\_ionized.xsc"];
set psfFile [glob "MtDARPin\_\_T3Q\_*\_\_solvated\_ionized.psf"];
set pdbFile [glob "MtDARPin\_\_T3Q\_*\_\_solvated\_ionized.pdb"];

#theta directory names that start with 0 will be interpreted as octals
proc forceInteger { x } {
    set count [scan $x %d%s n rest]
    if { $count <= 0 || ( $count == 2 && ![string is space $rest] ) } {
        return -code error "not an integer: \"$x\""
    }
    return $n
}

set fields [split $pdbFile "_"];
set z [lindex $fields 3];
set theta [lindex $fields 4];
#converts octal to integer
set theta [forceInteger $theta]
global z theta;

structure              		$psfFile
coordinates        			$pdbFile
outputName        			"o/MtDARPin\_\_T3Q\_$z\_$theta"

set temperature    300

# Continuing a job from the restart files
if {1} {
#binCoordinates     	HHC36_POPC-POPG_CHARMM36_equilibrate.coor
#binVelocities       	HHC36_POPC-POPG_CHARMM36_equilibrate.vel
extendedSystem	   	    $xscFile
}

firsttimestep      0


#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm	    on
parameters          ../par_all36_carb.prm
parameters          ../par_all36m_prot.prm
parameters			../par_all36_lipid.prm
parameters			../toppar_water_ions.mod.str
parameters          ../par_all36_na.prm
parameters          ../lig.prm
parameters          ../toppar_all36_na_nad_ppi.str
parameters          ../par_all36_cgenff.prm

# NOTE: Do not set the initial velocity temperature if you
# have also specified a .vel restart file!
temperature         $temperature


# Periodic Boundary conditions
# NOTE: Do not set the periodic cell basis if you have also
# specified an .xsc restart file!
if {0} {
cellBasisVector1     0.0    0.0		0.0
cellBasisVector2     0.0   	0.0   	0.0
cellBasisVector3     0.0	0.0		0.0
cellOrigin           0.0    0.0     0.0
}
wrapWater           on
wrapAll             on


# Force-Field Parameters
exclude             scaled1-4
1-4scaling          1.0
cutoff              12.
switching           on
switchdist          10.
pairlistdist        13.5
#margin 3.0


# Integrator Parameters
timestep            2.0  ;# 2fs/step
rigidBonds          all  ;# needed for 2fs steps
nonbondedFreq       1
fullElectFrequency  2
stepspercycle       10


#PME (for full-system periodic electrostatics)
if {1} {
PME                 	yes
PMEGridSpacing 		1.0
#PMEGridSizeX        24
#PMEGridSizeY        24
#PMEGridSizeZ        24
}


# Constant Temperature Control
if {1} {
langevin            on    ;# do langevin dynamics
langevinDamping     5     ;# damping coefficient (gamma) of 5/ps
langevinTemp        $temperature
langevinHydrogen    no    ;# don't couple langevin bath to hydrogens
}

# Constant Pressure Control (variable volume)
if {1} {
useGroupPressure      yes ;# needed for 2fs steps
useFlexibleCell       no  ;# no for water box, yes for membrane
useConstantArea       no  ;# no for water box, yes for membrane

langevinPiston        on
langevinPistonTarget  1.01325 ;#  in bar -> 1 atm
langevinPistonPeriod  100.
langevinPistonDecay   50.
langevinPistonTemp    $temperature
}


restartfreq         5000     ;# 5000steps = every 10ps
dcdfreq             5000
xstFreq             5000
outputEnergies      5000
outputPressure      5000


# Fixed restraints
if {1} {
constraints				on
consexp					2
conskfile				$pdbFile
consref					$pdbFile
conskcol				O
}




############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

if {1}	{
tclForces on

tclForcesScript {

proc sugar {first last} {
    set atoms {}
    for {set a $first} {$a <= $last} {incr a} {
        lappend  atoms $a;
        }
    return $atoms;
    }

set tcl_precision 17;
set PI  3.1415926535897931;

set dt 0.000002;
set v 20; # drift velocity of the nanobody waist in Å/ns.
# set v 2; # slow drift velocity of the ligand in Å/ns.
set tau  50100;  # Number of Time-steps for
set k 1 ; #Spring constant for peptide molecule in Kcal/(mol.Å^2)
set r0 25;
set n_first 2388;     # First atom of the sugar molecule
set n_last 2447;      # Last atom of the sugar molecule

set sug [sugar $n_first $n_last];
set T3Q [addgroup $sug];


############################################################################################################
############################################################################################################
############################################################################################################
proc calcforces {} {


set tcl_precision 17;

global v  PI  dt  r0  k  tau	z theta  T3Q;
loadcoords coor;
set t [getstep];	# i starts at 0, not 1, presumably as set by the firsttimestep variable


if {$t < $tau} { set radial $r0} else {
    set radial [expr $r0 - $v * ($t - $tau + 1) * $dt];}
set x [expr $radial*cos($theta*$PI/180)];
set y [expr $radial*sin($theta*$PI/180)];

set r "$x $y $z";
#set forces [open "forces\_$theta\_$phi.txt" "a"] ;
#puts $forces "[expr $i - $tau]		$r0		$coor($nanobody_waist)		$f";
#close $forces ;
# Finding and applying the forces, starting by finding the acceleration according to Verlet's algorithm
# Note: [vecsub $a $b] gives a-b
# set f  [ vecscale [ vecsub $r  $coor($T3Q) ] $k ];
set f  [ vecscale [ vecsub $r  $coor($T3Q) ] $k ];
addforce $T3Q $f;

}
#clearconfig;
}

}
#############################################################
## EXECUTION SCRIPT                                        ##
#############################################################

# Minimization
if {1} 		{
minimize            100
#reinitvels          $temperature
}

#run 50000000 ;# 100 ns
# run 550000 ;# 1.1 ns
# run 500000; #1 ns
run 445000; #0.85
#250000; #0.5 ns
