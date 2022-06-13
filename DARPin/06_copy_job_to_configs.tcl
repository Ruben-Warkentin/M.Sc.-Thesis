# Azimuthal angles around the z-axis, for a total of 14 segments
set thetas "015 045 075 105 135 165 195 225 255 285 315 345";
# z values at the center of each patch, from z = -22 to z = 27 Angstrom
set zz "-24.24 -17.31 -10.39 -3.46 03.46 10.39 17.31 24.24";

for {set i 0} {$i < [llength $zz]} { incr i }  {
	for {set j 0} {$j < [llength $thetas]} {incr j}	{
		## First, theta
		set theta [lindex $thetas $j]
		set z [lindex $zz $i];

        set job "configurations/-3.46\_015/job_SMD.sh";
        # set  "configurations/$z\_$theta";
        set destination "configurations/$z\_$theta";
        file copy -force $job $destination
		}
}
