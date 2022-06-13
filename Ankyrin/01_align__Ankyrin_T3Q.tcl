tkcon font "Courier bold" 12
set projection orthographic
light 0 on
light 1 on
light 2 on
light 3 on
axes location origin
color Display Background white

# Ankyrin structure taken from the last frame of an equilibration run, after a total of #### ns of equilibration at 300K and atmospheric pressure. 
# Water molecules and ions are then removed and the protein is moved so that
# its center of mass is situated at the origin.
# We start from that structure here:
mol new 4XCZ.psf
mol addfile 4XCZ.pdb

proc cartToCyl {r} {
    set pi 3.141592653589793
    set x [lindex $r 0]
    set y [lindex $r 1]
    set z [lindex $r 2]
    set rr [expr {sqrt($x*$x + $y*$y)}]
    if $x>0 { if $y>0 { set theta [expr {atan ($y/$x)*180/$pi}]} else {set theta [expr { 360 + atan ($y/$x)*180/$pi}]}}
    if $x<0 {set theta [expr {180 + atan ($y/$x)*180/$pi}]}
    if $x==0 {if $y>0 {set theta 90} else {set theta 270}}
    return "$rr $theta $z"
}

# Delete the default representation
# Documentation for the mol command: https://www.ks.uiuc.edu/Research/vmd/current/ug/node140.html
mol delrep 0 top
#mol selection {protein}
mol representation NewCartoon 0.3 30 4.1 0

# Create a base representation, showing all of the protein
mol selection {protein}
mol addrep top
mol modcolor 0 top Chain
mol modmaterial 0 top HardPlastic

# Add another representation to highlight the backbone in a different color:
mol color ResType
mol addrep top
mol modmaterial 1 top HardPlastic
mol modselect 1 top "protein and {resid 269 to 291 or resid 307 to 327 or resid 340 to 364 or resid 377 to 395}"

# Making two separate selections, one for the whole protein and one for its backbone:
set p [atomselect top "all"]
set bb [atomselect top "resid 269 to 291 or resid 307 to 327 or resid 340 to 364 or resid 377 to 395"]
set chains [atomselect top "protein and not (resid 269 to 291 or resid 307 to 327 or resid 340 to 364 or resid 377 to 395)"]

# Centering the backbone:
set c [center_of_mass $bb]
$p moveby [vecscale $c -1.00]
# Finding the symmetry axes
set inertia [measure inertia $bb]
set v3 [ lindex [lindex $inertia 1] 2]
# Draw v3 (scaled up) to see if it represents the symmetry axis of the backbone:
#draw color blue
#draw arrow "0 0 0" [vecscale $v3 25]
# Calculating M1, the 4x4 matrix that transforms the v3 vector to align parallel to the x-axis:
set M1 [transvecinv $v3]
# Performing the actual move:
$p move $M1
# Calculating the 4x4 matrix that rotates any vector by -90 degrees around the y-axis,
# so that v3 aligns with the z-axis:
set M2 [transaxis y -90]
# Performing the actual move:
$p move $M2
# Calculating the updated symmetry axes:
set inertia [measure inertia $bb]
set v3 [ lindex [lindex $inertia 1] 2]
# Centering the rotated molecule again, using the backbone coodinates:
set c [center_of_mass $p]
set z_c [lindex $c 2]
set c "0 0 $z_c"
$p moveby [vecscale $c -1.00]

# Measuring the spread of the chains across different dimensions:
# Get a list of chain residues:
set res_list [$chains get resid]
# Remove duplicates:
# set res_list [lsort -unique $res_list]
# set res_list [lsort $res_list]
puts $res_list

set res_r {}
set res_theta {}
set res_z   {}
for {set i 0} {$i < [llength $res_list]} {incr i} {
    set r [center_of_mass [atomselect top "protein and resid [lindex $res_list $i]"]]
    set r [cartToCyl $r]
    #puts $r
    lappend res_r [lindex $r 0]
    lappend res_theta [lindex $r 1]
    lappend res_z [lindex $r 2]
    }

set min_r [tcl::mathfunc::min {*}$res_r]
set max_r [tcl::mathfunc::max {*}$res_r]

set min_theta [tcl::mathfunc::min {*}$res_theta]
set max_theta [tcl::mathfunc::max {*}$res_theta]

set min_z [tcl::mathfunc::min {*}$res_z]
set max_z [tcl::mathfunc::max {*}$res_z]

puts "For the protein:"
puts "min r: $min_r   max r: $max_r"
puts "min θ: $min_theta   max θ: $max_theta"
puts "min z: $min_z   max z: $max_z"

## By Ruben
#For the side chains:
# min r: 5.711567715833719   max r: 20.047670333843275
# min θ: 204.44437299970906   max θ: 349.6710958001014
# min z: -14.69862247709215   max z: 24.07576117737416
# the radius determines the end point of the simulation, or the inner radius to
# which the ligand will be pushed to 

# finding the z-extend of the molcule
measure minmax $p
# but the lower z-extent of the whole molecule is about -26.4 Å, and the upper extent is about 27.8 Å.
# Let us assume we want to approch the backbone down to r = 5.7 Å
# Initially, we can bring each carbohydrate molecule from r = 25 Å inward, to have a total trip of 19.3 Å per simulation.
# We can easily go to even smaller r, but that likely introduces too many forced H-bonds that blur the main signal for our probing.

# Final destination of the carbohydrate molecules is on this cylinder:
draw material Transparent
draw color pink
draw cylinder {0 0 -26.5}  {0 0 27.9} radius 5.7 resolution 300
draw color blue
# And they each begin from the surface of this outer cylinder:
draw cylinder {0 0 -26.5}  {0 0 27.9} radius 25 resolution 300

# To be conservative, let us require that at r = 5.7*2 = 11.4 Å, we should have "full coverage", i.e., enough many carbohydrate probes
# that fully mosaic-cover the cylinder with r = 11.4 Å.
# Here is that cylinder:
draw color green
draw cylinder {0 0 -26.5}  {0 0 27.9} radius 11.4 resolution 300
# The surface area to cover is thus: A = 2πrh = 2π * 11.4 * (27.9+26.5) = 3896.6 Å^2
# Given that radius of gyration of T3Q is about 4.98 Å, each SLeX molecule can safely
# be considered a square mosaic of side 4.98*sqrt(2) = 7.04 ≈ 7 Å, havinge a 
# surface area of 49 Å^2. The number of simulation boxes is thus 3896.5/49 ≈ 79 boxes.
# The circumferance of the base of the cylinder at r = 11.4 Å is about 72 Å, whereas 
# the height of the cylinder is about 43.9 Å. So we need ~1.32 as many "mosaics"
# around the cylinder, as that along its height: 1.32*n*n ≈ 79, so n ≈ 8. 
# This way, we devide the circumferance of the cylinder into 10 segments and
# the height of the cylinder into 8 segments. This gives us 80 "mosaics".
