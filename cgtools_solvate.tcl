# +-------+
# | CG-it |
# +-------+

# CG-it, a VMD package to simplify creating Coarse Grained SDK 
# topologies.

# Copyright (c) 2013 by Chris MacDermaid <chris.macdermaid@gmail.com>
# and Giacomo Fiorin <giacomo.fiorin@gmail.com>

#(Shinoda) Shinoda, DeVane, Klein, Mol Sim, 33, 27 (2007).
#(DeVane) Shinoda, DeVane, Klein, Soft Matter, 4, 2453-2462 (2008).

## Solvate a selection in a waterbox with user specified padding
proc ::CGit::solvate_pad {sel {pad 10} {r 3.0} {waterbox default}} {

    variable datadir
    set molid [$sel molid]
    set seltext [$sel text]

    ## Load the waterbox
    if {$waterbox == "default"} {
        set waterbox $datadir/pdb/waterbox-cg.eq303
    }

    if {[catch {mol new $waterbox\.psf type psf} wb] ||
	[catch {mol addfile $waterbox\.pdb type pdb} wb]} {
        set msg [vmdcon -err "Unable to open $waterbox: $wb"]
        return -code error $msg
    }

    ## Copy the selection (solute) to a new mol
    set newmol [::TopoTools::selections2mol $sel]

    ## move the solute to the origin
    set sel2 [atomselect top "all"]
    set com [measure center $sel2 weight mass]
    $sel2 moveby [vecinvert $com]
    $sel2 set segname A79

    ## Dimensions of the solute + pad
    lassign [measure minmax $sel2] min max
    set min [vecadd $min [vecinvert [list $pad $pad $pad]]]
    set max [vecadd $max [list $pad $pad $pad]]
    set solute_dim [vecsub $max $min]

    ## center the waterbox and get its dimensions
    set wb_sel [atomselect $wb "all"]
    set com [measure center $sel2 weight mass]
    $wb_sel moveby [vecinvert $com]
    $wb_sel set segname B79

    # Dimensions of the solvent
    lassign [measure minmax $wb_sel] min max
    set solvent_dim [vecsub $max $min]

    ## Calculate fractional units of the solute w.r.t. the solvent box
    set frac {}
    foreach x $solute_dim y $solvent_dim {
        lappend frac [expr {ceil($x / $y)}]
    }

    ## Check if we need to replicate the box
    if {[vecsum $frac] > 3} {
        set wb_new [::TopoTools::replicatemol $wb {*}$frac]
        mol delete $wb
        set wb $wb_new
    }

    $sel2 delete
    $wb_sel delete

    ## Merge the solvent/solute
    set merged [::TopoTools::mergemols [list $newmol $wb]]

    ## Remove solvent near the solute
    set sel1 [atomselect $merged "same residue as (within $r of segname A79) and segname B79"]
    $sel1 set segname XX79
    $sel1 delete

    ## Trim the box to the minimum padded dimensions
    lassign $solute_dim a b c
    set aby2 [expr {$a / 2.0}]
    set bby2 [expr {$b / 2.0}]
    set cby2 [expr {$c / 2.0}]

    set dim_sel "x < [expr {-1 * $aby2}] or x > $aby2 or y < [expr {-1 * $bby2}] or y > $bby2 or z < [expr {-1 * $cby2}] or z > $cby2"
    set sel1 [atomselect $merged "($dim_sel) and segname B79"]
    $sel1 set segname XX79
    $sel1 delete

    set sel1 [atomselect $merged "not segname XX79"]
    set retval [::TopoTools::selections2mol $sel1]
    $sel1 delete

    ## Make sure we set the correct box dimensions
    molinfo $retval set {a b c} [list $a $b $c]

    ## Cleanup
    mol delete $newmol
    mol delete $wb
    mol delete $merged

    return $retval
}

## Fill the box based on the PBC dimensions from the solute molecule
proc ::CGit::solvate_box {sel {r 3.0} {waterbox default}} {

    variable datadir
    set molid [$sel molid]
    set seltext [$sel text]

    if {[$sel num] == 0} {
        set msg [vmdcon -err "Selection contains no atoms"]
        return -code error $msg
    }

    ## Load the waterbox
    if {$waterbox == "default"} {
        set waterbox $datadir/pdb/waterbox-cg.eq303
    }

    if {[catch {mol new $waterbox\.psf type psf} wb] ||
	[catch {mol addfile $waterbox\.pdb type pdb} wb]} {
        set msg [vmdcon -err "Unable to open $waterbox: $wb"]
        return -code error $msg
    }

    ## Copy the selection (solute) to a new mol
    set newmol [::TopoTools::selections2mol $sel]

    ## move the solute to the origin
    set sel2 [atomselect top "all"]
    set com [measure center $sel2 weight mass]
    $sel2 moveby [vecinvert $com]
    $sel2 set segname A79

    ## center the waterbox and get its dimensions
    set wb_sel [atomselect $wb "all"]
    set com [measure center $sel2 weight mass]
    $wb_sel moveby [vecinvert $com]
    $wb_sel set segname B79

    # Dimensions of the solvent
    lassign [measure minmax $wb_sel] min max
    set solvent_dim [vecsub $max $min]

    ##Get the pbc dimensions (to fill the box)
    lassign [molinfo $molid get {a b c}] a b c

    ## Calculate fractional units of the pbc-box w.r.t. the solvent box
    set frac {}
    foreach x [list $a $b $c] y $solvent_dim {
        lappend frac [expr {ceil($x / $y)}]
    }

    ## Check if we need to replicate the box
    if {[vecsum $frac] > 3} {
        set wb_new [::TopoTools::replicatemol $wb {*}$frac]
        mol delete $wb
        set wb $wb_new
    }

    $sel2 delete
    $wb_sel delete

    ## Merge the solvent/solute
    set merged [::TopoTools::mergemols [list $newmol $wb]]

    ## Remove solvent near the solute
    set sel1 [atomselect $merged "same residue as (within $r of segname A79) and segname B79"]
    $sel1 set segname XX79
    $sel1 delete

    ## Trim the waterbox to the pbc-box
    set aby2 [expr {$a / 2.0}]
    set bby2 [expr {$b / 2.0}]
    set cby2 [expr {$c / 2.0}]

    set dim_sel "x < [expr {-1 * $aby2}] or x > $aby2 or y < [expr {-1 * $bby2}] or y > $bby2 or z < [expr {-1 * $cby2}] or z > $cby2"
    set sel1 [atomselect $merged "($dim_sel) and segname B79"]
    $sel1 set segname XX79
    $sel1 delete

    set sel1 [atomselect $merged "not segname XX79"]
    set retval [::TopoTools::selections2mol $sel1]
    $sel1 delete

    ## Make sure we set the correct box dimensions
    molinfo $retval set {a b c} [list $a $b $c]

    ## Cleanup
    mol delete $newmol
    mol delete $wb
    mol delete $merged

    return $retval
}
