
def script_rg_rmsd(cord: str, top: str, traj: str, prefix_out:str, path: str, selection:str) -> str:
    script_rg_rmsd = """
    package require pbctools 

    mol load """ + cord.split(".")[-1] + " " + cord + " " + top.split(".")[-1] + " " + top + """ 
    mol addfile """ +  traj + """ type """ +  traj.split(".")[-1] + """waitfor all

    pbc unwrap -sel  '""" + selection + """'

    set outfile [open """ +  path + "/" + prefix_out + """.dat w];
    set nf [molinfo top get numframes]
    set sel_protein_cg [atomselect top '""" + selection + """']
    set frame0 [atomselect top '""" + selection + """' frame 0]
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

    return script_rg_rmsd


def script_rg_rmsd_big_traj(cord: str, top: str, traj: str, prefix_out:str, path: str, selection:str) -> str:
    script_rg_rmsd_big_traj = """
    package require pbctools 
    source bigdcd.tcl

    mol load """ + cord.split(".")[-1] + " " + cord + " " + top.split(".")[-1] + " " + top + """  

    set outfile [open """ +  path + "/" + prefix_out + """.dat w];
    set sel_protein_cg [atomselect top '""" + selection + """']
    set frame0 [atomselect top '""" + selection + """' frame 0]
    puts $outfile "RG\tRMSD"

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
        
    pbc unwrap -sel  '""" + selection + """'
    $sel_protein_cg move [measure fit $sel_protein_cg $frame0]
    puts $outfile "[gyr_radius $sel_protein_cg]\t[measure rmsd $sel_protein_cg $frame0]"
    }

    bigdcd rg_rmsd """ + traj + """
    bigdcd_wait

    close $outfile
    quit
    """
    return script_rg_rmsd_big_traj

