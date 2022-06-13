# Azimuthal angles around the z-axis, for a total of 14 segments
set thetas "015 045 075 105 135 165 195 225 255 285 315 345";
# z values at the center of each patch, from z = -22 to z = 27 Angstrom
set zz "-25.33 -18.09 -10.86 -3.62 03.62 10.86 18.09 25.33";

for {set i 0} {$i < [llength $zz]} { incr i }  {
	for {set j 0} {$j < [llength $thetas]} {incr j}	{
		## First, theta
		set theta [lindex $thetas $j]
		set z [lindex $zz $i];

        set SMD_run "configurations/-3.62\_015/SMD_run.tcl";
        # set  "configurations/$z\_$theta";
        set destination "configurations/$z\_$theta";
        # file copy -force configurations/-4.5\_0.00/SMD_run.tcl configurations/$z\_$theta
        file copy -force $SMD_run $destination
		}
}
