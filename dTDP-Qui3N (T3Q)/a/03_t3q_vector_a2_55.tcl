#this scripts creates a vector between two atoms and writes the unit vector over time 

# tkcon font "courier" 14
# display projection orthographic
# axes location origin

#loading dcd and pdb file from the simulation
mol new t3q_dry.dcd waitfor all
mol addfile t3q_x0.pdb

#defining numframes
set nf [molinfo top get numframes]
set outfile [open v_list_a2_55.txt w]
#puts $outfile "Frame# A2A_v"

	for {set i 0} {$i < $nf} {incr i} {
	animate goto $i
	#defining a vector between two atoms
	#the vector
	#create a vector between atoms index 0 and index 84
		#select the two atoms
		set a2 [atomselect top "index 2"]
		set a55 [atomselect top "index 55"]

	 	#get coordinates
		set a2_coor [lindex [$a2 get {x y z}] 0]
		set a55_coor [lindex [$a55 get {x y z}] 0]

		set v [vecsub $a55_coor $a2_coor]
		set lv [veclength $v]
		set uv [vecscale [expr 1/$lv] $v]
	puts $outfile "$uv"
	}
	close $outfile