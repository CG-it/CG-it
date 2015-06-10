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

# +----------------------------------------------------------------------+
# | Measure all bonds, angles, dihedrals, impropers and store them in    |
# | array arr based on types: arr(bond a-b), arr(angle a-b-c), arr(dihed |
# | a-b-c-d), arr(imprp a-b-c-d)                                         |
# +----------------------------------------------------------------------+

proc ::CGtools::traj_analyze {sel {start 0} {end -1} {flags {bonds angles}} } {

    variable sys

    ## Array to store the calculated distributions in
    catch {array unset parr *}
    array set parr {}

    set molid [$sel molid]
    set seltext [$sel text]
    set nframes [molinfo $molid get numframes]

    if {$end == -1} {set end [expr {$nframes - 1}]}

    foreach f $flags {
        switch -exact -- $f {
            bonds     {
                ## Measure all bonds
                topo -sel $sel retypebonds
                set bondlist [topo -sel $sel getbondlist type]
                vmdcon -info "measuring bonds"
                foreach b $bondlist {
                    lassign $b id1 id2 type
                    lappend parr([list bond $type])\
                        [measure bond [list $id1 $id2] molid $molid first $start last $end]
                }

            }

            angles    {
                ## Measure all angles
                topo -sel $sel guessangles
                set anglelist [topo -sel $sel getanglelist]
                vmdcon -info "measuring angles"
                foreach a $anglelist {
                    lassign $a type id1 id2 id3
                    lappend parr([list angle $type])\
                        [measure angle [list $id1 $id2 $id3] molid $molid first $start last $end]
                }
            }

            dihedrals {
                ## Measure all Dihedrals
                topo -sel $sel guessdihedrals
                set dihedrallist [topo -sel $sel getdihedrallist]
                vmdcon -info "measuring dihedrals"
                foreach d $dihedrallist {
                    lassign $d type id1 id2 id3 id4
                    lappend parr([list dihed $type])\
                        [measure dihed [list $id1 $id2 $id3 $id4] molid $molid first $start last $end]
                }
            }

            impropers {
                ## Measure all Impropers
                topo -sel $sel guessimpropers tolerance 45; # Tol should be large enough to get them all...
                set improperlist [topo -sel $sel getimproperlist]
                vmdcon -info "measuring impropers"
                foreach i $improperlist {
                    lassign $i type id1 id2 id3 id4
                    lappend parr([list imprp $type])\
                        [measure imprp [list $id1 $id2 $id3 $id4] molid $molid first $start last $end]
                }

            }
        }
    }

    ## Print out the data
    traj_analyze_print parr

    return $sys(OK)
}

proc ::CGtools::traj_analyze_print {arr} {

    variable sys

    ## Array to store the calculated distributions in
    upvar $arr parr

    ## Type of measurement
    foreach m {bond angle dihed imprp} {

        ## Measure and type
        foreach x [lsort -ascii -index 1 [array names parr "$m*"]] {
            lassign $x label type

            set fid [open [format "%s_%s.dat" $m $type] "w"]
            puts $fid [format "##SDK %s %s" $m $type]

            ## Each individual
            foreach y $parr($x) {

                ## Each Frame
                foreach z $y {
                    puts $fid [format "%11.6f" $z]
                }

            }

            close $fid

        }
    }

    return $sys(OK)
}

proc ::CGtools::traj_density {molid {start 0} {end -1} {stride 1}} {

    set nframes [molinfo $molid get numframes]
    if {$end == -1} {set end [expr {$nframes - 1}]}


    puts "Using [expr {int(($end - $start) / $stride) + 1}] frames to calculate density"

    ## Mass (g/mol)
    set sel [atomselect $molid "all"]
    set mass [vecsum [$sel get mass]]
    $sel delete

    ## Check for zero mass
    if {[expr {$mass == 0.0}]} {return 0.0}

    set V {}
    for {set i $start} {$i <= $end} {incr i $stride} {
        ##Volume (\A**3)
        lassign [molinfo $molid get {a b c}] a b c
        lappend V [expr {$a * $b * $c}]
    }

    ## Average the volume
    set V [vecmean $V]

    ## Return the density in g/mL
    return [expr {$mass / $V * 1.660539}]
}
