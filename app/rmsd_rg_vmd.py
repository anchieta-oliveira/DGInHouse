import os


def _tcl_quote(value: str) -> str:
    text = os.fspath(value)
    for character in ('\\', '"', '$', '[', ']'):
        text = text.replace(character, '\\' + character)
    return f'"{text}"'


def script_rg_rmsd(cord: str, top: str, traj: str, prefix_out: str, path: str, selection: str) -> str:
    cord_type = os.path.splitext(cord)[1].lstrip('.')
    top_type = os.path.splitext(top)[1].lstrip('.')
    traj_type = os.path.splitext(traj)[1].lstrip('.')
    return f'''
package require pbctools

mol load {cord_type} {_tcl_quote(cord)} {top_type} {_tcl_quote(top)}
mol addfile {_tcl_quote(traj)} type {traj_type} waitfor all

pbc unwrap -sel {_tcl_quote(selection)}

set outfile [open {_tcl_quote(os.path.join(path, prefix_out + '.dat'))} w]
set nf [molinfo top get numframes]
set sel_protein_cg [atomselect top {_tcl_quote(selection)}]
set frame0 [atomselect top {_tcl_quote(selection)} frame 0]
puts $outfile "RG\\tRMSD"

proc gyr_radius {{sel}} {{
    if {{[$sel num] <= 0}} {{
        error "gyr_radius: must have at least one atom in selection"
    }}
    set com [measure center $sel weight mass]
    set sum 0
    foreach coord [$sel get {{x y z}}] {{
        set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
    }}
    return [expr sqrt($sum / ([$sel num] + 0.0))]
}}

for {{set i 1}} {{$i < $nf}} {{incr i}} {{
    $sel_protein_cg frame $i
    $sel_protein_cg move [measure fit $sel_protein_cg $frame0]
    puts $outfile "[gyr_radius $sel_protein_cg]\\t[measure rmsd $sel_protein_cg $frame0]"
}}
close $outfile
quit
'''
