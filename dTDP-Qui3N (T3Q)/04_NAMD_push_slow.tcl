#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# t3q being pushed in one direction -- SMD

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################
set PSFFile 01_t3q_x0_ionized.psf;
set PDBFile 01_t3q_x0_ionized.pdb;

structure              		$PSFFile
coordinates        			$PDBFile
outputName        			"o_slow/t3q_push_slow_01"

set temperature    300

# Continuing a job from the restart files
if {0} {
#binCoordinates     	HHC36_POPC-POPG_CHARMM36_equilibrate.coor
#binVelocities       	HHC36_POPC-POPG_CHARMM36_equilibrate.vel
extendedSystem	   	    $XSCFile
}

firsttimestep      0


#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm	    on
parameters          par_all36_carb.prm
parameters          par_all36m_prot.prm
parameters			par_all36_lipid.prm
parameters			toppar_water_ions.mod.str
parameters          par_all36_na.prm
parameters          lig.prm
parameters          toppar_all36_na_nad_ppi.str
parameters          par_all36_cgenff.prm

# NOTE: Do not set the initial velocity temperature if you
# have also specified a .vel restart file!
temperature         $temperature


# Periodic Boundary conditions
# NOTE: Do not set the periodic cell basis if you have also
# specified an .xsc restart file!
if {1} {
cellBasisVector1 50.457 0 0 
cellBasisVector2 0 40.226 0 
cellBasisVector3 0 0 43.304 
cellOrigin -19.4085 -0.8410 -0.4900 
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
# PMEGridSizeX        41
# PMEGridSizeY        40
# PMEGridSizeZ        40
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
dcdfreq             5000     ;# 500 steps = every 1ps
xstFreq             5000
outputEnergies      5000
outputPressure      5000


# Fixed restraints
if {0} {
constraints				on
consexp					2
conskfile				$PDBFile
consref					$PDBFile
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
        lappend atoms $a;
        }
    return $atoms;    
}

set tcl_precision 17;
set PI  3.1415926535897931;

set dt 0.000002;
# set v 20; # drift velocity of the nanobody waist in Å/ns.
set v 2; # drift velocity of the nanobody waist in Å/ns.
set tau 100;  # Number of Time-steps for minimization + equilibration
set k 1 ; #Spring constant for peptide molecule in Kcal/(mol.Å^2)
set x0 0;
set n_first 1;     # First atom of the sugar molecule
set n_last 60;      # Last atom of the sugar molecule

set sug [sugar $n_first $n_last];
set t3q [addgroup $sug];

############################################################################################################
############################################################################################################
############################################################################################################
proc calcforces {} {


set tcl_precision 17;

global v  PI  dt  x0  k  tau  t3q;
loadcoords coor;
set t [getstep];	# i starts at 0, not 1, presumably as set by the firsttimestep variable


if {$t < $tau} { set x $x0} else {
    set x [expr $x0 + $v * ($t - $tau + 1) * $dt];}
set y [expr 0];
set z [expr 0];

set r "$x $y $z";
#set forces [open "forces\_$theta\_$phi.txt" "a"] ;
#puts $forces "[expr $i - $tau]		$r0		$coor($nanobody_waist)		$f";
#close $forces ;
# Finding and applying the forces, starting by finding the acceleration according to Verlet's algorithm
# Note: [vecsub $a $b] gives a-b
# set f  [ vecscale [ vecsub $r  $coor($SLeX) ] $k ];
set f  [ vecscale [ vecsub $r  $coor($t3q) ] $k ];
addforce $t3q $f;

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

run 50000000 ;#100ns