#script for recording hbonds between protein and ligand
# hbonds -sel1 [atomselect top "chain H"] -sel2 [atomselect top "chain N and resid 99 to 109"] -writefile yes -type all -outfile $Hbonds_file_CDR3   -detailout $details_file_CDR3

package require hbonds

# Azimuthal angles around the z-axis, for a total of 12 segments
set thetas "015	045	075	105	135	165	195	225	255	285	315	345";
set thetas_dcd "15	45	75	105	135	165	195	225	255	285	315	345";
# z values at the center of each patch, from z = -27.7 to z = 27.7 Angstrom
set zz "-24.2375 -17.3125 -10.3875 -3.4625 3.4625 10.3875 17.3125 24.2375";
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
        set z [format "%05.2f" $z]

        set psfFile "configurations/$z\_$theta/2XEE\_\_T3Q\_$z\_$theta\_\_solvated\_ionized";
        set dcdFile "configurations/$z\_$theta/o/2XEE\_\_T3Q\_$z\_$theta_dcd";
        mol new $psfFile.psf
        mol addfile $dcdFile.dcd waitfor all
        hbonds -sel1 [atomselect top "segid 2XEE"] -sel2 [atomselect top "segid T3Q"] -outdir configurations/$z\_$theta/  -writefile yes -type all
        mol delete all;
        }
    }
