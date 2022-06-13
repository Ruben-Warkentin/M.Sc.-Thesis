########################################################################
########################################################################
package require solvate
package require autoionize

# Azimuthal angles around the z-axis, for a total of 12 segments
set thetas "15 45 75 105 135 165 195 225 255 285 315 345";
# z values at the center of each patch, from z = -27.7 to z = 27.7 Angstrom
set zz "-25.33 -18.09 -10.86 -3.62 3.62 10.86 18.09 25.33";

for {set i 0} {$i < [llength $zz]} { incr i }  {
	for {set j 0} {$j < [llength $thetas]} {incr j}	{
		## First, theta
		set theta [lindex $thetas $j]
		set z [lindex $zz $i];
		set theta [format "%03.0f" $theta]
		set z [format "%05.2f" $z]

		set psfFile "configurations/$z\_$theta/MtDARPin\_\_T3Q\_$z\_$theta.psf";
		set pdbFile "configurations/$z\_$theta/MtDARPin\_\_T3Q\_$z\_$theta.pdb";
		solvate $psfFile $pdbFile -t 20 -o "configurations/$z\_$theta/MtDARPin\_\_T3Q\_$z\_$theta\_\_solvated";
		file delete -force $psfFile $pdbFile

		set psfFile "configurations/$z\_$theta/MtDARPin\_\_T3Q\_$z\_$theta\_\_solvated.psf";
		set pdbFile "configurations/$z\_$theta/MtDARPin\_\_T3Q\_$z\_$theta\_\_solvated.pdb";
		set logFile "configurations/$z\_$theta/MtDARPin\_\_T3Q\_$z\_$theta\_\_solvated.log";
		autoionize -psf $psfFile -pdb $pdbFile -sc 0.15 -o "configurations/$z\_$theta/MtDARPin\_\_T3Q\_$z\_$theta\_\_solvated\_ionized"
		mol delete all;
		file delete -force $psfFile $pdbFile $logFile


		file mkdir "configurations/$z\_$theta/o";
		set xscFile "configurations/$z\_$theta/MtDARPin\_\_T3Q\_$z\_$theta\_\_solvated\_ionized.xsc";
		set psfFile "configurations/$z\_$theta/MtDARPin\_\_T3Q\_$z\_$theta\_\_solvated\_ionized.psf";
		set pdbFile "configurations/$z\_$theta/MtDARPin\_\_T3Q\_$z\_$theta\_\_solvated\_ionized.pdb";
		mol new $psfFile;
		mol addfile $pdbFile;
		set a [atomselect top all]
		$a set beta 0;
		$a set occupancy 0;
		set bb [atomselect top "protein and resid 1 to 25 or resid 38 to 58 or resid 71 to 91 or resid 104 to 124 or resid 137 to 157"];
		$bb set occupancy 1;
		$a writepdb $pdbFile;
		$a writepsf $psfFile;
		mol delete all;

		mol new $psfFile;
		mol addfile $pdbFile;
		set a [atomselect top all];
		set mm [measure minmax $a];
		set d [vecsub [lindex $mm 1] [lindex $mm 0]];
		set c [measure center $a];
		set xsc [open $xscFile "w"];
		puts $xsc "# NAMD extended system configuration restart file";
		puts $xsc "#\$LABELS step a\_x a\_y a\_z b\_x b\_y b\_z c\_x c\_y c\_z o\_x o\_y o\_z";
		puts $xsc "0 [lindex $d 0] 0 0 0 [lindex $d 1] 0 0 0 [lindex $d 2] $c";
		close $xsc;
		mol delete all;
		}
}
