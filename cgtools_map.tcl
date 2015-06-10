#!/usr/bin/tclsh

# +----------+
# | CG-Tools |
# +----------+

# CGtools, a VMD package to simplify creating coarse grained SDK
# topologies.

# Copyright (c) 2013 by Chris MacDermaid <chris.macdermaid@gmail.com>
# and Giacomo Fiorin <giacomo.fiorin@gmail.com>

#(Shinoda) Shinoda, DeVane, Klein, Mol Sim, 33, 27 (2007).
#(DeVane) Shinoda, DeVane, Klein, Soft Matter, 4, 2453-2462 (2008).

proc CGtools::map {args} {

    set sel [lindex $args 0]
    set molid [$sel molid]
    set seltext [$sel text]
    set seltext "(same residue as ($seltext))"
    set frame [molinfo $molid get frame]

    set args [lrange $args 1 end]
    set nargs [llength $args]
    set newargs ""

    ## Default Arguments
    set first $frame; set last $frame; set stride 1;

    ## Make a local selection that contains entire residues,
    ## we can only CG entire residues for now.
    set localsel [atomselect $molid $seltext]

    ## Check for flags and remove them from arguments
    for {set i 0} {$i < $nargs} {incr i} {
        set arg [lindex $args $i]
        if {[string match -?* $arg]} {
            set val [lindex $args [expr $i+1]]
            switch -- $arg {
                -first  {set first  $val}
                -last   {set last   $val}
                -stride {set stride $val}
                -all    {set first 0; set last\
                             [molinfo $molid get numframes]}
                -now     {set first $frame; set last $frame}
                default {cgCon -info "Map: Unknown option $arg"; return}
            }
        } else {
            lappend newargs $arg
        }
    }

    ## Change to the first frame
    molinfo $molid set frame $first

    ## Change the occupancy to a unique beadID
    ## This beadID is used when mapping to distinguish
    ## separate beads.
    setBeadID $localsel

    ## Get the order inwhich we should traverse VMD's
    ## coordinate array such that atoms belonging to the
    ## same bead are grouped together
    setOrder $localsel

    ## We coarsegrain one third of the water in the system.
    ## Identify those molecules using a unique identifier.
    ## This unfortunately destroys any dynamical information
    ## in the system w.r.t. water.
    $sel set user4 1; # > 0 == include atom in mapping
    set sel2 [atomselect $molid "$seltext and water and name \"O.*\""]
    set nwater [$sel2 num]
    if {$nwater > 0} {
        set sel3 [atomselect $molid "water"]; $sel3 set user4 -1; $sel3 delete;
        set user4 {}
        set residues [$sel2 get residue]
        foreach {r1 r2 r3} $residues {
            lappend user4 1 -1 -1
        }
        $sel2 set user4 [lrange $user4 0 $nwater-1]
    }
    $sel2 delete

    set props {}
    set keys {residue resid resname chain segname segid occupancy beta}

    ## Get unique residues defined in the selection, currently only entire
    ## residues can be coarse-grained
    set sel2 [atomselect $molid "$seltext and user4 > 0"]
    set residues [lsort -unique -increasing -integer -index 0\
                      [$sel2 get $keys]]

    ## Create the beads for each residue
    foreach r $residues {
        lassign $r {*}$keys
        foreach bead [getMap $resname -keys {type name charge mass} -transpose] {
            lappend props [concat $bead $chain $resid $resname\
                               $segname $segid $occupancy $beta]
        }
    }

    set keys {type name charge mass chain\
                  resid resname segname segid occupancy beta}

    ## Make a new "empty" mol
    set molname CGmol-$molid
    set newmol -1
    set ntotal [llength $props]
    if {[catch {mol new atoms $ntotal} newmol]} {
        cgCon -error "map: could not create new molecule: $molid"
        cgCon -error $newmol
        return -1
    }

    mol rename $newmol $molname

    ## Set the properties
    set sel3 [atomselect $newmol "all"]
    $sel3 set $keys $props
    $sel3 delete

    ## Update PBC
    molinfo $newmol set {a b c} [molinfo $molid get {a b c}]

    ## Update internal VMD structures
    mol reanalyze $newmol

    ## Map (Will automatically use the c-routine if the library is loaded)
    cgmap -sel $sel2 -append $newmol -first $first -last $last -stride $stride

    ## Cleanup
    $sel2 delete
    $localsel delete

    ## Return the created MOLID
    return $newmol
}

## Assigns a unique BeadID to atoms
## using the user field.
proc CGtools::setBeadID {sel} {

    set molid   [$sel molid]
    set seltext [$sel text]

    ## Get the residues specified in the selection
    set resnames [lsort -unique -ascii [$sel get resname]]

    ## Clear the user field
    $sel set user -1.0

    foreach r $resnames {
        set sel2 [atomselect $molid "($seltext) and resname $r"]

        ## Get a list of beads for this residue
        set beads [getMap $r -keys {map}]
        set nbeads($r) [llength $beads]

        ## Create lookup hash
        set id 0
        foreach b $beads {
            foreach a $b {
                set lookup($r,$a) $id
            }
            incr id
        }
        $sel2 delete
    }

    ## Multiply the unique residue id by
    ## the lookup index
    set ids {}
    foreach name [$sel get name]\
        resid [$sel get residue]\
        resn  [$sel get resname] {

            if {[catch {set idx $lookup($resn,$name)}]} {
                cgCon -err "Unknown atom $name in residue type $resn"
                return -code error
            }

            lappend ids [expr {int($resid * $nbeads($resn) + $idx)}]
        }
    $sel set user $ids
}

proc CGtools::setOrder {sel} {

    set molid   [$sel molid]
    set seltext [$sel text]
    set seltext "(($seltext) and not water)"

    ## Get the residues specified in the selection
    set resnames [lsort -unique -ascii [$sel get resname]]

    ## Clear the user fields
    $sel set user2 [$sel get index]
    $sel set user3 [$sel get mass]

    ## This is not the most efficient or elegant, but it's
    ## fast enough and it gets the job done.
    ## Sort by bead order supplied in the map, then pass that
    ## to the mapper in the "user" fields.
    foreach r $resnames {
        set sel2 [atomselect $molid "($seltext) and resname $r"]
        set segid [$sel2 get segid]
        set beads [getMap $r -keys {map}]
        $sel2 set {segid user4 user user2 user3}\
            [lsort -command [list aacompare [join $beads]]\
                 [$sel2 get {name residue user index mass}]]
        $sel2 set segid $segid
        $sel2 delete
    }
}

## Compares residue first, then
## by particular atom order specified
## as a list in "l"
proc CGtools::aacompare {l a b} {

    lassign $a name_a residue_a
    lassign $b name_b residue_b

    if {$residue_a < $residue_b} {
        return -1

    } elseif {$residue_a > $residue_b} {
        return 1

    } else {

        set i [lsearch -ascii -exact $l $name_a]
        set j [lsearch -ascii -exact $l $name_b]

        return [expr {$i < $j ? -1 : $i > $j ? 1 : 0}]
    }
}

## Map the coordinates
## The best spead you can get with this routine
## is about 20 fps.
proc CGtools::cgmap {args} {

    set nargs [llength $args]
    set newargs ""

    ## Default Arguments
    set first 0; set last 0; set stride 1;
    set molid [molinfo top]; set appendmol -1; set sel ""

    ## Check for flags and remove them from arguments
    for {set i 0} {$i < $nargs} {incr i} {
        set arg [lindex $args $i]
        if {[string match -?* $arg]} {
            set val [lindex $args [expr $i+1]]
            switch -- $arg {
                -molid  {set molid     $val}
                -append {set appendmol $val}
                -first  {set first     $val}
                -last   {set last      $val}
                -stride {set stride    $val}
                -sel    {set sel $val; \
                             set molid [$sel molid]}
                -all    {set first 0; set last\
                             [molinfo $molid get numframes]}
                -now     {set first $frame; set last $frame}
                default {cgCon -info "Map: Unknown option $arg"; return}
            }
        } else {
            lappend newargs $arg
        }
    }

    if {$molid < 0 || $appendmol < 0} {
        cgCon -err "Must specify molids to work with"
        return -code error
    }

    set localsel 0
    if {$sel == ""} {
        set sel [atomselect $molid "all"]
        set localsel 1
    }

    set nframes [molinfo $molid get numframes]
    set natoms  [$sel num]
    set sel_all [atomselect $molid "all"]
    set sel2 [atomselect $appendmol "all"]

    ## Output to screen
    set print [expr {($last - $first) / 10}]
    if {$print < 10} {set print 10}
    if {$print > 10} {set print 100}

    ## BeadIDs and array access order
    set bid [$sel get user]; set order [$sel get user2]
    lappend bid -1

    ## Mass is invariant over frames
    set m [$sel_all get mass];

    for {set frame $first}\
        {$frame <= $last && $frame < $nframes} {incr frame $stride} {

            if {[expr {$frame % $print}] == 0} {
                cgCon -info "Mapping frame $frame"
            }

            ## Update the frame/selections
            molinfo $molid set frame $frame
            $sel_all update; $sel update

            set x [$sel_all get x]; set y [$sel_all get y]; set z [$sel_all get z]

            set com {};
            set w 0; set cx 0; set cy 0; set cz 0;

            for {set i 0} {$i < $natoms} {incr i} {
                set idx [expr int([lindex $order $i])]
                set b1 [lindex $bid $i]
                set b2 [lindex $bid $i+1]

                set xx [lindex $x $idx]
                set yy [lindex $y $idx]
                set zz [lindex $z $idx]
                set mm [lindex $m $idx]

                if {$b1 != $b2} {
                    set w  [expr {$w + $mm}]
                    set cx [expr {$cx + ($mm * $xx)}]
                    set cy [expr {$cy + ($mm * $yy)}]
                    set cz [expr {$cz + ($mm * $zz)}]
                    lappend com [vecscale\
                                     [list $cx $cy $cz]\
                                     [expr {1.0/$w}]]
                    set w 0; set cx 0; set cy 0; set cz 0;
                } else {
                    set w  [expr {$w + $mm}]
                    set cx [expr {$cx + ($mm * $xx)}]
                    set cy [expr {$cy + ($mm * $yy)}]
                    set cz [expr {$cz + ($mm * $zz)}]
                }
            }

            ## Make a new frame and set coordinates
            animate dup $appendmol; $sel2 update
            $sel2 set {x y z} $com

            ## Update the PBCs
            molinfo $appendmol set {a b c}\
                [molinfo $molid get {a b c}]
        }

    if {$localsel} {$sel delete}
    $sel2 delete
    $sel_all delete
}
