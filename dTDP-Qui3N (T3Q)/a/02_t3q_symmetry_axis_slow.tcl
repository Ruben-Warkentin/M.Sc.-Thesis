# tkcon font "Courier bold" 14
display projection orthographic
axes location origin

#loading dcd and pdb file from the simulation
mol new t3q_dry_slow.dcd waitfor all
mol addfile t3q_x0.pdb
#setting nf variable to the # of frames
set nf [molinfo top get numframes]
#set an outfile, with write access
set outfile_v1 [open 02_slow_symmetry_axis_v1.txt w]
set outfile_v2 [open 02_slow_symmetry_axis_v2.txt w]
set outfile_v3 [open 02_slow_symmetry_axis_v3.txt w]
#write into the outfile, there are headings
puts $outfile_v1 "symmetry_v1_slow x_y_z"
puts $outfile_v2 "symmetry_v2_slow x_y_z"
puts $outfile_v3 "symmetry_v3_slow x_y_z"
#create a loop
#this will contain the script that calcualtes the symmetry axis
#at each frame, and then store the unit vector that represents then
#symmetry axis
for {set i 0} {$i < $nf} {incr i} {
    #advance frame
    animate goto $i
    #select the molecule, SLeX in this case
    set s [atomselect top "not water and not ions"]
    #add frame xxx after "" if you want to select a specific frame

    #axis of symmetry
    set inertia [measure inertia $s]
    set v1 [ lindex [lindex $inertia 1] 0]
    set v2 [ lindex [lindex $inertia 1] 1]
    set v3 [ lindex [lindex $inertia 1] 2]

    ##drawing the new inertia vector
    # draw color purple
    # draw line "0 0 0" [vecscale $v3 25]

    #this stores the vector in the outfile txt file
    puts $outfile_v1 [lindex $v1]
    puts $outfile_v2 [lindex $v2]
    puts $outfile_v3 [lindex $v3]
}
close $outfile_v1
close $outfile_v2
close $outfile_v3