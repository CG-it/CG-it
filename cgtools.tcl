#!/usr/bin/tclsh

# +-------+
# | CG-it |
# +-------+

# CG-it, a VMD package to simplify creating coarse grained SDK
# topologies.

# Copyright (c) 2013 by Chris MacDermaid <chris.macdermaid@gmail.com>
# and Giacomo Fiorin <giacomo.fiorin@gmail.com>

#(Shinoda) Shinoda, DeVane, Klein, Mol Sim, 33, 27 (2007).
#(DeVane) Shinoda, DeVane, Klein, Soft Matter, 4, 2453-2462 (2008).

## None of this will work without topotools.
package require topotools

## This package makes extensive use of some Tcl 8.5 features.
package require Tcl 8.5

namespace eval ::CGit:: {

    namespace export CGit

    variable version 1.1

    ## Valid cgtools commands

    variable validcmd {readtop readparam reanalyzemol setbonds
        setangles setdihedrals reset clear writelammpsdata writelammpsparam
        writesdktopo writetop setbondparam setangleparam setpairparam
        deletebondparam deleteangleparam deletepairparam getbondparam
        getangleparam getpairparam readmap map maptrajectory maptraj
        printmap setparam deleteparam getparam printparam printmaptable
        printmap solvate help commands topo2map}

    ## Parameter Files loaded by default
    variable param_files {param_sdk.dat}

    ## Parameter array
    #variable params
    #array set params {}
    #array unset params *

    #Parameter Dictionary
    variable params ""

    ## Store interface variables
    variable sys
    array unset sys

    ## Atom Mappings and Topologies
    variable map ""

    ## Test map internal consistency on loading?
    variable maptest 0

    ## Show the title and citation info?
    variable showtitle 0

    ## Locations of additional data files
    variable datadir $env(CGTOOLSDIR)

    ## Make some utility commands accessable
    namespace export map map_trajectory

    ## Return error codes
    set sys(OK) 0
    set sys(ERROR) -1
}

# +----------------+
# | Global command |
# +----------------+

proc CGit { args } {
    eval ::CGit::CGit $args
}

proc CGit::title {} {

    variable version

    cgCon -info "CG-it Version $version    "
    cgCon -info "Loading from: [info script] "
    cgCon -info "Please Cite:                "

}

proc CGit::usage {} {

    set groups(1) "common flags:"
    set commands(1) {\
                         { {-molid}       {<num>|top}         {molecule id (default: 'top')}}\
                         { {-sel}         {<selection>}       {atom selection function or text (default: 'all')}}\
                     }

    set groups(2) ""
    set commands(2) {\
                         { {help} {} {prints this message}}\
                     }

    set groups(3) "Topologies and Parameters:"
    set commands(3) {\
                         { {readtop}       {<filename>}                 {read topology file}}\
                         { {readparam}     {<filename>}                 {read parameter file}}\
                         { {writetop}      {[<filename>]}               {write legacy SDK topology for specified selection}}\
                         { {topo2map}      {<topfile> [<mapfile>]}      {convert legacy SDK topology to map}}\
                         { {printparam}    {}                           {print to screen all currently loaded parameters}}\
                         { {clear}         {<maps> <params> <both>}     {clear the all maps, params or both}}\
                         { {reset}         {}                           {clear loaded parameters and topologies, reload defaults}}\
                         { {} {} {}}\
                         { {(set|delete)bondparam}  {<type1> <type2>}            {re/set bond parameters for specified atom types}}\
                         { {}              {<k> <r_bond>}               {}}\
                         { {(set|delete)angleparam} {<type1> <type2> <type3>}    {re/set angle parameters for specified atom types}}\
                         { {}              {<k> <theta>}                {}}\
                         { {(set|delete)pairparam}  {<type1> <type2> <lj_type>}  {re/set pair parameters for specified atom types}}\
                         { {}              {<epsilon> <sigma>}          {}}\
                         { {} {} {}}\
                         { {get(bond|angle|pair)param} {}               {print specified parameters. Wildcard '*' returns }}\
                         { {}              {<type1> <type2> [<type3> <type4>]} {all types: e.g. 'cg getbondparam PH *'}}\
                     }

    set groups(4) "CG mappings:"
    set commands(4) {\
                         { {readmap}       {<filename>}             {read an SDK map file}}\
                         { {map}           {}                       {perfrom a CG mapping}}\
                         { {maptraj}       {[<start> <end> <skip>]} {perfrom a CG mapping of a trajectory}}\
                         { {printmap}      {<resname>}              {print out mapping information for specified}}\
                         { {}              {}                       {residue name (default 'all')}}\
                     }

    set groups(5) "Molecule Preparation:"
    set commands(5) {\
                         {  {set(bonds|angles|dihedrals)} {}          {Apply bonding/angle connectivities specified from}}\
                         {  {}               {}                        {topologies/maps to selected atoms.}}\
                         {  {reanalyzemol}   {}                        {Apply [<mass> <charge> <index> <type> <all> <none>] properties}}\
                         {  {}               {}                        {from topologies to selected atoms and retype}}\
                         {  {}               {}                        {[<bonds> <angles> <dihedrals> <impropers>] (default 'all bonds angles')}}\
                         {  {solvate}        {[<distance> <filename>]} {Solvate selected atoms in a solvent box according}}\
                         {  {}               {}                        {to pbc cell size with a minimum <distance> (default '3')}}\
                         {  {}               {}                        {angstroms from the solute. <filename> specifies a pdb file}}\
                         {  {}               {}                        {for a orthorhombic solvent box.}}\

    }

    set groups(6) "LAMMPS:"
    set commands(6) {\
                         { {writelammpsdata}  {<filename> <atomstyle>} {Write out a data file suitable for use with lammps}}\
                         { {writelammpsparam} {<filename>}             {Write out a lammps SDK 'parameter' file containing}}\
                         { {}                 {<cpu> <gpu> <both>}     {forcefield parameters required by the system. Compatible}}\
                         { {}                 {}                       {options for simulating with <cpu> <gpu> or <both> are provided}}\
                     }

    cgCon -info "usage: cg <command> \[args...\] <flags>"
    cgCon -info "Note: by default, cgtools works on all atoms in the top molecule"
    cgCon -info ""

    ## Print out each group
    foreach g {2 1 3 4 5 6} {
        cgCon -info $groups($g)
        foreach cmd $commands($g) {cgCon -info [format "  %-15s %-25s %-25s" {*}$cmd]}
        cgCon -info ""
    }
}

# +---------------------------------------------------+
# | Parse arguments, set defaults, execute directives |
# +---------------------------------------------------+

proc CGit::CGit { args } {

    variable validcmd

    set molid -1
    set seltxt all
    set selmol -1
    set localsel 1

    set cmd {}
    set sel {}

    set newargs {}

    for {set i 0} {$i < [llength $args]} {incr i} {
        set arg [lindex $args $i]

        if {[string match -?* $arg]} {

            set val [lindex $args [expr $i+1]]

            switch -- $arg {
                -molid {
                    if {[catch {molinfo $val get name} res]} {
                        cgCon -err "Invalid -molid argument '$val': $res"
                        return
                    }
                    set molid $val
                    if {[string equal $molid "top"]} {
                        set molid [molinfo top]
                    }
                    incr i
                }

                -sel {
                    # check if the argument to -sel is a valid atomselect command
                    if {([info commands $val] != "") && ([string equal -length 10 $val atomselect])} {
                        set localsel 0
                        set selmol [$val molid]
                        set sel $val
                    } else {
                        set localsel 1
                        set seltxt $val
                    }
                    incr i
                }

                -- break

                default {
                    #cgCon -info "default: $arg"
                    lappend newargs $arg
                }
            }
        } else {
            lappend newargs $arg
        }
    }

    if {$molid < 0} {
        set molid $selmol
    }

    if {$molid < 0} {
        set molid [molinfo top]
    }

    set retval ""
    if {[llength $newargs] > 0} {
        set cmd [lindex $newargs 0]
        set newargs [lrange $newargs 1 end]
    } else {
        set newargs {}
        set cmd help
    }

    if {[lsearch -exact -ascii $validcmd $cmd] < 0} {
        cgCon -err "Unknown command '$cmd'"
        usage
        return
    }

    if { ![string equal $cmd help] && $molid >= 0 } {
        if {($selmol >= 0) && ($selmol != $molid)} {
            cgCon -err "Molid from selection '$selmol' does not match -molid argument '$molid'"
            return
        }

        if {$localsel} {
            # need to create a selection
            if {[catch {atomselect $molid $seltxt} sel]} {
                cgCon -err "Problem with atom selection using '$seltxt': $sel"
                return
            }
        }
    }

    switch -- $cmd {

        readtop {

            if {[llength $newargs] < 1} {
                cgCon -err "cgtools 'readtop' requires a file to read from"
                usage
                return
            }

            ## Convert legacy top to map on the fly
            readTopo {*}$newargs
        }

        topo2map {
            if {[llength $newargs] < 1} {
                usage
                return
            }

            set out [lindex $newargs 1]
            if {$out == ""} {set out stdout}

            topo2map [lindex $newargs 0] $out 0
        }

        writetop -
        writesdktopo {
            if {$sel == "" || $molid < 0} {return [usage]}
            writesdktopo $sel {*}$newargs
        }

        readparam {
            if {[llength $newargs] < 1} {
                cgCon -err "cgtools 'readparam' requires a file to read from"
                usage
                return
            }

            #read_param $newargs
            readParam {*}$newargs
        }

        writelammpsdata {

            if {[llength $newargs] < 1} {
                cgCon -err "cgtools 'writelammpsdata' requires a file to write to as argument"
                usage
                return
            }

            writelammpsdata $molid {*}$newargs
        }

        writelammpsparam {

            if {[llength $newargs] < 1} {
                cgCon -err "cgtools 'writelammpsparam' requires a file to write to as argument"
                usage
                return
            }

            writelammpsparam $molid {*}$newargs
        }

        printmap -
        printmaptable {map_stats {*}$newargs}

        mapread -
        readmap {
            if {[llength $newargs] < 0} {
                cgCon -err "cgtools 'readmap' requires a file to read from as argument"
                usage
                return
            }

            readTopo{*}$newargs
        }

        map -
        maptraj -
        maptrajectory {
            if {$sel == "" || $molid < 0} {return [usage]}
            set retval [map $sel {*}$newargs]
        }

        setbond -
        setbonds  {
            if {$sel == "" || $molid < 0} {return [usage]}
            set_bonds $sel
        }

        setangle -
        setangles {
            if {$sel == "" || $molid < 0} {return [usage]}
            set_angles $sel
        }

        setdihedral -
        setdihedrals {
            if {$sel == "" || $molid < 0} {return [usage]}
            set_dihedrals $sel
        }

        reanalyzemol {
            if {$sel == "" || $molid < 0} {return [usage]}
            reanalyze_mol $sel {*}$newargs
        }

        reset {reset}

        clear {
            clear {*}$newargs
        }

        solvate {
            if {$sel == "" || $molid < 0} {return [usage]}
            set retval [solvate_box $sel {*}$newargs]
        }

        setbondparam  {set_param bond  {*}$newargs}
        setangleparam {set_param angle {*}$newargs}
        setpairparam  {set_param pair  {*}$newargs}
        setparam      {set_param {*}$newargs}

        deletebondparam  {delete_param bond  {*}$newargs}
        deleteangleparam {delete_param angle {*}$newargs}
        deletepairparam  {delete_param pair  {*}$newargs}
        deleteparam   {delete_param {*}$newargs}

        getbondparam  {set retval [get_param bond {*}$newargs]}
        getangleparam {set retval [get_param angle {*}$newargs]}
        getpairparam  {set retval [get_param pair {*}$newargs]}
        getparam      {set retval [get_param {*}$newargs]}
        printparam    {print_param {*}$newargs}

        commands {return [join $validcmd]}

        help -
        default {usage}
    }

    if {$localsel} {
        catch {$sel delete}
    }

    return $retval
}

# +---------+
# | Helpers |
# +---------+

## A wrapper around cgCon
proc CGit::cgCon {flag str} {
    switch -- $flag {
        "-error" -
        "-err"  { vmdcon -err  "CG> $str" }
        "-warn" { vmdcon -warn "CG> $str" }
        "-info" { vmdcon -info "CG> $str" }
        default { vmdcon $str}
    }
}

# +---------+
# | Startup |
# +---------+

interp alias {} cg {} ::CGit::CGit
package provide cg $::CGit::version

if {$::CGit::showtitle} {::CGit::title}

## Load other cgtools files
#source [file join $env(CGTOOLSDIR) maps common.tcl]
source [file join $env(CGTOOLSDIR) cgtools_map.tcl]
source [file join $env(CGTOOLSDIR) cgtools_partopo.tcl]
source [file join $env(CGTOOLSDIR) cgtools_traj.tcl]
source [file join $env(CGTOOLSDIR) cgtools_lammps.tcl]
source [file join $env(CGTOOLSDIR) cgtools_solvate.tcl]
source [file join $env(CGTOOLSDIR) cgtools_legacy.tcl]
source [file join $env(CGTOOLSDIR) json.tcl]

## C-based routines
catch {load [file join $env(CGTOOLSDIR) libcgmap.so]}

## Load the default mappings from maps in $env(CGTOOLSDIR)/maps
##::CGit::map_defaults
