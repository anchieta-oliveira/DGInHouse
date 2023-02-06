SCRIPT_RG_RMSD= """
package require pbctools


mol load gro ../3oty_cg_ion.gro psf ../3oty_cg_ion.psf 
mol addfile 3oty_cg_md.xtc type xtc waitfor all


pbc unwrap -sel "resid 1 to 99"


set outfile [open rg_cg_and_rmsd.dat w];
set nf [molinfo top get numframes]
set sel_protein_cg [atomselect top "resid 1 to 99"]
set frame0 [atomselect top "resid 1 to 99" frame 0]
puts $outfile "RG\tRMSD"


proc gyr_radius {sel} {
  # make sure this is a proper selection and has atoms
  if {[$sel num] <= 0} {
    error "gyr_radius: must have at least one atom in selection"
  }
  # gyration is sqrt( sum((r(i) - r(center_of_mass))^2) / N)
  set com [measure center $sel weight mass]
  set sum 0
  foreach coord [$sel get {x y z}] {
    set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
  }
  return [expr sqrt($sum / ([$sel num] + 0.0))]
}


# rg and RMSD calculation loop
for {set i 1 } {$i < $nf } { incr i } {
    $sel_protein_cg frame $i
    $sel_protein_cg move [measure fit $sel_protein_cg $frame0]
    puts $outfile "[gyr_radius $sel_protein_cg]\t[measure rmsd $sel_protein_cg $frame0]"
}
close $outfile
quit """

SCRIPT_RG_RMSD_BIG_TRAJ= """
package require pbctools 
source bigdcd.tcl

mol load gro ../3OUC_cg_ion.gro psf ../3OUC_cg_ion.psf 

set outfile [open rg_cg_and_rmsd_com_lig_teste.dat w];
set sel_protein_cg [atomselect top "resid 1 to 99"]
set frame0 [atomselect top "resid 1 to 99" frame 0]
puts $outfile "RG\tRMSD"


# rg and RMSD calculation loop

proc rg_rmsd { frame } {
  global frame0 sel_protein_cg outfile
	  proc gyr_radius {sel} {
	  # make sure this is a proper selection and has atoms
	  if {[$sel num] <= 0} {
	    error "gyr_radius: must have at least one atom in selection"
	  }
	  # gyration is sqrt( sum((r(i) - r(center_of_mass))^2) / N)
	  set com [measure center $sel weight mass]
	  set sum 0
	  foreach coord [$sel get {x y z}] {
	    set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
	  }
	  return [expr sqrt($sum / ([$sel num] + 0.0))]
	}
	
   pbc unwrap -sel "resid 1 to 99"
   $sel_protein_cg move [measure fit $sel_protein_cg $frame0]
   puts $outfile "[gyr_radius $sel_protein_cg]\t[measure rmsd $sel_protein_cg $frame0]"
}

bigdcd rg_rmsd  3OUC_cg_md.xtc
bigdcd_wait

close $outfile
quit
"""
