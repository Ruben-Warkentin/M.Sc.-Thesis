#script for recording hbonds between protein and ligand
# hbonds -sel1 [atomselect top "chain H"] -sel2 [atomselect top "chain N and resid 99 to 109"] -writefile yes -type all -outfile $Hbonds_file_CDR3   -detailout $details_file_CDR3

package require hbonds

# Azimuthal angles around the z-axis, for a total of 12 segments
set thetas "018 054 090 126 162 198 234 270 306 342";
set thetas_dcd "18 54 90 126 162 198 234 270 306 342";
# z values at the center of each patch, from z = -27.7 to z = 27.7 Angstrom
set zz "-23.80 -17.00 -10.20 -3.40 03.40 10.20 17.00 23.80";
######################################################
#testing script with a single theta and z value

# 		## First, theta
# 		set theta 0.00;
# 		set z -4.5;
#
#         set psfFile "configurations/$z\_$theta/Loop\_DARPIN\_\_SLeX\_$z\_$theta\_\_solvated_ionized";
#         set dcdFile "configurations/$z\_$theta/o/Loop\_DARPIN\_\_SLeX\_$z\_$theta";
#         mol new $psfFile.psf
#         mol addfile $dcdFile.dcd waitfor all
#         hbonds -sel1 [atomselect top "protein"] -sel2 [atomselect top "segname SLeX"] -outdir configurations/$z\_$theta/  -writefile yes -type all
#         mol delete all;

######################################################
for {set i 0} {$i < [llength $zz]} { incr i }  {
	for {set j 0} {$j < [llength $thetas]} {incr j}	{
		## First, theta
		set theta [lindex $thetas $j]
        set theta_dcd [lindex $thetas_dcd $j]
		set z [lindex $zz $i];

        set psfFile "configurations/$z\_$theta/4XCZ\_\_T3Q\_$z\_$theta\_\_solvated\_ionized";
        set dcdFile "configurations/$z\_$theta/o/4XCZ\_\_T3Q\_$z\_$theta_dcd";
        mol new $psfFile.psf
        mol addfile $dcdFile.dcd waitfor all
        hbonds -sel1 [atomselect top "segid 4XCZ"] -sel2 [atomselect top "segid T3Q"] -outdir configurations/$z\_$theta/  -writefile yes -type all
        mol delete all;
        }
    }
