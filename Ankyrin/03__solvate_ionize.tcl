########################################################################
########################################################################
package require solvate
package require autoionize

# Azimuthal angles around the z-axis, for a total of 12 segments
# set thetas "18 54 90 126 162 198 234 270 306 342";
# z values at the center of each patch, from z = -27.7 to z = 27.7 Angstrom
# set zz "-23.80 -17.00 -10.20 -3.40 3.40 10.20 17.00 23.80";

# testing to remove solvating error that causes atom overlap
set thetas "18 54 90 126 162 198 234 270 306 342";
set zz "-23.80 17.00 -10.20 -3.40 3.40 10.20 17.00 23.80";


for {set i 0} {$i < [llength $zz]} { incr i }  {
	for {set j 0} {$j < [llength $thetas]} {incr j}	{
		## First, theta
		set theta [lindex $thetas $j]
		set z [lindex $zz $i];
		set theta [format "%03.0f" $theta]
		set z [format "%05.2f" $z]

		set psfFile "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta.psf";
		set pdbFile "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta.pdb";
		solvate $psfFile $pdbFile -t 20 -o "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta\_\_solvated";
		file delete -force $psfFile $pdbFile

		set psfFile "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta\_\_solvated.psf";
		set pdbFile "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta\_\_solvated.pdb";
		set logFile "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta\_\_solvated.log";
		autoionize -psf $psfFile -pdb $pdbFile -sc 0.15 -o "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta\_\_solvated\_ionized"
		mol delete all;
		file delete -force $psfFile $pdbFile $logFile


		file mkdir "configurations/$z\_$theta/o";
		set xscFile "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta\_\_solvated\_ionized.xsc";
		set psfFile "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta\_\_solvated\_ionized.psf";
		set pdbFile "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta\_\_solvated\_ionized.pdb";
		mol new $psfFile;
		mol addfile $pdbFile;
		set a [atomselect top all]
		$a set beta 0;
		$a set occupancy 0;
		set bb [atomselect top "protein and resid 269 to 291 or resid 307 to 327 or resid 340 to 364 or resid 377 to 395"];
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
