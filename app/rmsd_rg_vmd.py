
def script_rg_rmsd(cord: str, top: str, traj: str, prefix_out:str, path: str, selection:str) -> str:
    script_rg_rmsd = """
    package require pbctools 

    mol load """ + cord.split(".")[-1] + " " + cord + " " + top.split(".")[-1] + " " + top + """ 
    mol addfile """ +  traj + """ type """ +  traj.split(".")[-1] + """ waitfor all

    pbc unwrap -sel  " """ + selection + """ "

    set outfile [open """ +  path + "/" + prefix_out + """.dat w];
    set nf [molinfo top get numframes]
    set sel_protein_cg [atomselect top " """ + selection + """ "]
    set frame0 [atomselect top " """ + selection + """ " frame 0]
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
    set sel_protein_cg [atomselect top " """ + selection + """ "]
    set frame0 [atomselect top " """ + selection + """ " frame 0]
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
        
    pbc unwrap -sel  " """ + selection + """ "
    $sel_protein_cg move [measure fit $sel_protein_cg $frame0]
    puts $outfile "[gyr_radius $sel_protein_cg]\t[measure rmsd $sel_protein_cg $frame0]"
    }

    bigdcd rg_rmsd """ + traj + """
    bigdcd_wait

    close $outfile
    quit
    """
    return script_rg_rmsd_big_traj

def script_bigdcd() -> str:
    bigdcd = """
    # BigDCD-v2
    # Justin Gullingsrud
    # updates by Axel Kohlmeyer
    # vmd@ks.uiuc.edu
    ################################################################################
    # Purpose: Use this script to analyze one or more trajectory files that 
    # don't fit into memory.  The script will arrage for your analysis function 
    # to be called each time a frame is loaded, then delete the frame from memory.
    # The analysis script must accept one argument; BigDCD will keep track of how
    # many timesteps have been loaded and call your script with that number.
    #
    # How to include this function: either source the script directory, or 
    # (better) place the script in one of the directories in your auto_path 
    # variable and include "package require bigdcd" in your script.
    #
    # New in version 2:
    #
    # The bigdcd command accepts an (optional) argument that defines
    # the format type of the trajectory file. Its default is "auto", i.e.
    # let VMD guess it from the file name.
    #
    # The function bigdcd_wait can be called to have a script wait until
    # all frames have been processed. This is most useful for batch processing
    # of analysis script. By construction bigdcd executes in the background 
    # and will return the control to VMD after the last trajectory file has 
    # been scheduled for reading. Thus analysis scripts will terminate prematurely.
    ################################################################################

    # Example 1: 
    # This computes the center of mass for each frame in the DCD file.
    #
    # proc mycenter { frame } {
    #   global all
    #   puts "$frame: [measure center $all weight mass]"
    # }
    # set $mol [mol new alanin.psf type psf waitfor all]
    # set all [atomselect $mol all]
    # $all global
    # bigdcd mycenter auto alanin.dcd 
    #

    # Example 2:
    # This computes the RMS distance between each frame in
    # a sequence of xyz files and a reference pdb file.  
    # this example even works in batch mode.
    #
    # proc myrmsd { frame } {
    #   global ref sel all
    #   $all move [measure fit $sel $ref]
    #   puts "$frame: [measure rmsd $sel $ref]"
    # }
    # set mol [mol new protein.psf type psf waitfor all]
    # set all [atomselect $mol all]
    # set ref [atomselect $mol "name CA" frame 0]
    # set sel [atomselect $mol "name CA"]
    # mol addfile protein.pdb waitfor all
    # bigdcd myrmsd xyz eq01.xyz eq02.xyz eq03.xyz
    # bigdcd_wait
    # quit

    proc bigdcd { script type args } {
        global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame bigdcd_running
    
        set bigdcd_running 1
        set bigdcd_frame 0
        set bigdcd_firstframe [molinfo top get numframes]
        set bigdcd_proc $script

        # backwards "compatibility". type flag is omitted.
        if {[file exists $type]} { 
            set args [linsert $args 0 $type] 
            set type auto
        }
    
        uplevel #0 trace variable vmd_frame w bigdcd_callback
        foreach dcd $args {
            if { $type == "auto" } {
                mol addfile $dcd waitfor 0
            } else {
                mol addfile $dcd type $type waitfor 0
            }
        }
        after idle bigdcd_wait
    }

    proc bigdcd_callback { tracedvar mol op } {
        global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame
        set msg {}
    
        # If we're out of frames, we're also done 
        # AK: (can this happen at all these days???). XXX
        set thisframe $vmd_frame($mol)
        if { $thisframe < $bigdcd_firstframe } {
            puts "end of frames"
            bigdcd_done
            return
        }
    
        incr bigdcd_frame
        if { [catch {uplevel #0 $bigdcd_proc $bigdcd_frame} msg] } { 
            puts stderr "bigdcd aborting at frame $bigdcd_frame\n$msg"
            bigdcd_done
            return
        }
        animate delete beg $thisframe end $thisframe $mol
        return $msg
    }

    proc bigdcd_done { } {
        global bigdcd_running
        
        if {$bigdcd_running > 0} then {
            uplevel #0 trace vdelete vmd_frame w bigdcd_callback
            puts "bigdcd_done"
            set bigdcd_running 0
        }
    }

    proc bigdcd_wait { } {
        global bigdcd_running bigdcd_frame
        while {$bigdcd_running > 0} {
            global bigdcd_oldframe
            set bigdcd_oldframe $bigdcd_frame
            # run global processing hooks (including loading of scheduled frames)
            display update ui
            # if we have read a new frame during then the two should be different.
            if { $bigdcd_oldframe == $bigdcd_frame } {bigdcd_done}
        }
    }
    """
    return bigdcd