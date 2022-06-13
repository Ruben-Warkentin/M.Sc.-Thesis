This text file contains info on the rotational autocorrelation analysis
of t3q

1. CatDCD for getting only the sugar atoms
 a) create in index list of atoms 0-59
 b) tcl script:
 set filename "t3q_0_59"
 set file [open $filename.ind w]
 foreach indices [$sel get index] {
 puts $file "$indices"
 }
 close $file

c) creating a dcd with select atoms 
catdcd -o output.dcd -i index.ind trajectory.dcd

1) "fast" simulation
catdcd -o a/t3q_dry.dcd -i a/t3q_0_59.ind o_x0/t3q_push_01.dcd

2) "slow" simulation
catdcd -o a/t3q_dry_slow.dcd -i a/t3q_0_59.ind o_slow/t3q_push_slow_01.dcd
