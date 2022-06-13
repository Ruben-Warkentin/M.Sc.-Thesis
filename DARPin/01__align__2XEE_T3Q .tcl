tkcon font "Courier bold" 12
set projection orthographic
light 0 on
light 1 on
light 2 on
light 3 on
axes location origin
color Display Background white

# 2XEE DARPin structure taken from the last frame of a 2021_02_09__DARPin--2XEE__02__51.86_ns, after a total of 106.18 ns of equilibration. 
# Water molecules and ions are then removed and the DARPin molecule is moved so that
# its center of mass is situated at the origin.
# We start from that structure here:
mol new 2XEE.psf
mol addfile 2XEE.pdb

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
mol modselect 1 top "protein and {resid 13 to 37 or resid 50 to 70 or resid 83 to 103 or resid 116 to 135 or resid 149 to 169}"

# Making two separate selections, one for the whole protein and one for its backbone:
set p [atomselect top "all"]
set bb [atomselect top "resid 13 to 37 or resid 50 to 70 or resid 83 to 103 or resid 116 to 135 or resid 149 to 169"]
set chains [atomselect top "protein and not (resid 13 to 37 or resid 50 to 70 or resid 83 to 103 or resid 116 to 135 or resid 149 to 169)"]

set c [center_of_mass $bb]
$p moveby [vecscale $c -1.00]
# Finding the symmetry axes
set inertia [measure inertia $bb]
set v3 [ lindex [lindex $inertia 1] 2]
# Draw v3 (scaled up) to see if it represents the symmetry axis of the backbone:
#draw color blue
#draw arrow "0 0 0" [vecscale $v3 25]
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
# Centering the rotated molecule again, using the backbone coodinates:
set c [center_of_mass $bb]
$p moveby [vecscale $c -1.00]
draw color purple
draw arrow "0 0 0" [vecscale $v3 35]

# Measuring the spread of the chains across different dimensions:
# Get a list of chain residues:
set res_list [$chains get resid]
# Remove duplicates:
set res_list [lsort -unique $res_list]
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

puts "For the side chains:"
puts "min r: $min_r   max r: $max_r"
puts "min θ: $min_theta   max θ: $max_theta"
puts "min z: $min_z   max z: $max_z"

# Output:
# For the side chains:
# min r: 6.392893555623027   max r: 19.45463312879902
# min θ: 0.9315940815926248   max θ: 356.59499330341475
# min z: -19.902887929639014   max z: 21.057607234840713

# but the lower z-extent of the whole molecule is about -27.2 Å, and the upper extent is about 28.2 Å.
# Let us assume we want to approch the backbone down to r = 6.4 Å
# Initially, we can bring each carbohydrate molecule from r = 25 Å inward, to have a total trip of 18.6 Å per simulation.
# We can easily go to even smaller r, but that likely introduces too many forced H-bonds that blur the main signal for our probing.

# Final destination of the carbohydrate molecules is on this cylinder:
draw color pink
draw cylinder {0 0 -27.2}  {0 0 28.2} radius 6.4 resolution 300
draw color blue
# And they each begin from the surface of this outer cylinder:
draw cylinder {0 0 -27.2}  {0 0 28.2} radius 25 resolution 300

# To be conservative, let us require that at r = 6.4*2 = 12.8 Å, we should have "full coverage", i.e., enough many carbohydrate probes
# that fully mosaic-cover the cylinder with r = 12.8 Å.
# Here is that cylinder:
draw color green
draw cylinder {0 0 -27.2}  {0 0 28.2} radius 12.8 resolution 300
# The surface area to cover is thus: A = 2πrh = 2π * 12.8 * (27.2+28.2) = 4455 Å^2
# Given that radius of gyration of T3Q is about 4.98 Å, each SLeX molecule can safely
# be considered a square mosaic of side 4.98*sqrt(2) = 7.04 ≈ 7 Å, havinge a 
# surface area of 49 Å^2. The number of simulation boxes is thus 4455/49 ≈ 91 boxes.
# The circumferance of the base of the cylinder at r = 12.8 Å is about 80.4 Å, whereas 
# the height of the cylinder is about 55.4 Å. So we need ~1.45 as many "mosaics"
# around the cylinder, as that along its height: 1.45*n*n ≈ 91, so n ≈ 8. 
# This way, we devide the circumferance of the cylinder into 11 segments and
# the height of the cylinder into 8 segments.
# Again, we can be even more conservative and use 12 segments along the circumferance, to get
# a total of 96 segments (rather than 88).
