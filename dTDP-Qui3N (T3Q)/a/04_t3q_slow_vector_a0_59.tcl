#this scripts creates a vector between two atoms and writes the unit vector over time 

# tkcon font "courier" 14
# display projection orthographic
# axes location origin

#loading dcd and pdb file from the simulation
mol new t3q_dry_slow.dcd waitfor all
mol addfile t3q_x0.pdb

#defining numframes
set nf [molinfo top get numframes]
set outfile [open v_list_a0_59_slow.txt w]
#puts $outfile "Frame# A2A_v"

	for {set i 0} {$i < $nf} {incr i} {
	animate goto $i
	#defining a vector between two atoms
	#the vector
	#create a vector between atoms index 0 and index 84
		#select the two atoms
		set a0 [atomselect top "index 0"]
		set a59 [atomselect top "index 59"]

	 	#get coordinates
		set a0_coor [lindex [$a0 get {x y z}] 0]
		set a59_coor [lindex [$a59 get {x y z}] 0]

		set v [vecsub $a59_coor $a0_coor]
		set lv [veclength $v]
		set uv [vecscale [expr 1/$lv] $v]
	puts $outfile "$uv"
	}
	close $outfile
