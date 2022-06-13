# Azimuthal angles around the z-axis, for a total of 14 segments
set thetas "018 054 090 126 162 198 234 270 306 342";
# z values at the center of each patch, from z = -22 to z = 27 Angstrom
set zz "-23.80 -17.00 -10.20 -3.40 03.40 10.20 17.00 23.80";

for {set i 0} {$i < [llength $zz]} { incr i }  {
	for {set j 0} {$j < [llength $thetas]} {incr j}	{
		## First, theta
		set theta [lindex $thetas $j]
		set z [lindex $zz $i];

        set job "configurations/-3.40\_018/job_SMD.sh";
        # set  "configurations/$z\_$theta";
        set destination "configurations/$z\_$theta";
        # file copy -force configurations/-4.5\_0.00/SMD_run.tcl configurations/$z\_$theta
        file copy -force $job $destination
		}
}
