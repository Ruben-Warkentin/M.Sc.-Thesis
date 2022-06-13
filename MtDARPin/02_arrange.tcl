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

mol new MtDARPin_T3Q.psf;
mol addfile MtDARPin_T3Q.pdb;
set pi 3.14159265359;

# Azimuthal angles around the z-axis, for a total of 12 segments
set thetas "15 45 75 105 135 165 195 225 255 285 315 345";
# z values at the center of each patch, from z = -27.7 to z = 27.7 Angstrom
set zz "-25.33 -18.09 -10.86 -3.62 3.62 10.86 18.09 25.33";

set r 25; # Initial radius (distance from the z-zxis)

set all [atomselect top all]
$all set beta 0
$all set occupancy 0
set mtdarpin [atomselect top "segid MT"]
set T3Q [atomselect top "segid T3Q"]

set bb [atomselect top "protein and resid 1 to 25 or resid 38 to 58 or resid 71 to 91 or resid 104 to 124 or resid 137 to 157"];
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
		$all writepdb "configurations/$z\_$theta/MtDARPin__T3Q_$z\_$theta.pdb";
		$all writepsf "configurations/$z\_$theta/MtDARPin__T3Q_$z\_$theta.psf";
		}
}
