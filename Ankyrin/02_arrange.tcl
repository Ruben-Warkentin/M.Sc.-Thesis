proc reset {m}	{
	# resets the location and orientation of the molecule,
	# such that its center of mass is at the origin and its p1 aligns parallel to the z-axis.
	
	# Move to the origin
	set cm [center_of_mass $m];
	$m moveby [vecscale -1.000 $cm];
	
	set inertia [measure  inertia $m];
	set p1 [ lindex [lindex $inertia 1] 0];

	# Align with the x-axis
	set M [transvecinv $p1];
	$m move $M;
		
	# Now align against the z-axis
	set M [transabout "0 1 0" 90];
	$m move $M;

	# Move to the origin
	set cm [center_of_mass $m];
	$m moveby [vecscale -1.000 $cm];
}
########################################################################
########################################################################

mol new 4XCZ_T3Q.psf;
mol addfile 4XCZ_T3Q.pdb;
set pi 3.14159265359;

# Azimuthal angles around the z-axis, for a total of 12 segments
set thetas "18 54 90 126 162 198 234 270 306 342";
# z values at the center of each patch, from z = -27.7 to z = 27.7 Angstrom
set zz "-23.80 -17.00 -10.20 -3.40 3.40 10.20 17.00 23.80";

set r 25; # Initial radius (distance from the z-zxis)

set all [atomselect top all]
$all set beta 0
$all set occupancy 0
set ankyrin [atomselect top "segid 4XCZ"]
set T3Q [atomselect top "segid T3Q"]

set bb [atomselect top "protein and resid 269 to 291 or resid 307 to 327 or resid 340 to 364 or resid 377 to 395"];
set p [atomselect top "protein"]
$bb set occupancy 1;


set c [center_of_mass $p]
$p moveby [vecscale $c -1.00]
# Finding the symmetry axes
set inertia [measure inertia $bb]
set v3 [ lindex [lindex $inertia 1] 2]
# Calculating M, the 4x4 matrix that transforms the v3 vector to align parallel to the x-axis:
set M1 [transvecinv $v3]
# Performing the actual move:
$p move $M1
# Calculating the 4x4 matrix that rotates any vector by -90 degrees around the y-axis
set M2 [transaxis y -90]
# Performing the actual move:
$p move $M2
# Calculating the updated symmetry axes:
set inertia [measure inertia $bb]
set v3 [ lindex [lindex $inertia 1] 2]
# Centering the rotated molecule one last time, just along the z-axis, based on the spread of the whole protein:
set c [center_of_mass $p]
set z_c [lindex $c 2]
set c "0 0 $z_c"
$p moveby [vecscale $c -1.00]

reset $T3Q;
$T3Q set beta 1;

for {set i 0} {$i < [llength $zz]} { incr i }  {
	for {set j 0} {$j < [llength $thetas]} {incr j}	{
		reset $T3Q;
		## First, z
		set z [lindex $zz $i];
		## Then theta
		set theta [lindex $thetas $j];
		set x [expr $r*cos($theta*$pi/180)]
		set y [expr $r*sin($theta*$pi/180)]

        $T3Q moveby "$x $y $z"
		
		set theta [format "%03.0f" $theta]
		set z [format "%05.2f" $z]
		
		file mkdir "configurations/$z\_$theta";
		$all writepdb "configurations/$z\_$theta/4XCZ__T3Q_$z\_$theta.pdb";
		$all writepsf "configurations/$z\_$theta/4XCZ__T3Q_$z\_$theta.psf";
		}
}
