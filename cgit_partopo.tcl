# +-------+
# | CG-it |
# +-------+

# CG-it, a VMD package to simplify creating Coarse Grained SDK
# topologies.

# Copyright (c) 2013 by Chris MacDermaid <chris.macdermaid@gmail.com>
# and Giacomo Fiorin <giacomo.fiorin@gmail.com>

#(Shinoda) Shinoda, DeVane, Klein, Mol Sim, 33, 27 (2007).
#(DeVane) Shinoda, DeVane, Klein, Soft Matter, 4, 2453-2462 (2008).

## Read param wrapper based on provided or guessed filetype
proc CGit::readParam {filename {type guess}} {

    cgCon -info "Reading parameters from $filename"

    if {![string compare $type "guess"]} {
        set type [file extension $filename]
    }

    switch -exact $type {
        ".json" -
        "json"  {set retval [readParamJSON $filename]}
        ".dat"  -
        "dat"   {set retval [readParamOLD $filename]}
        default {cgCon -err "Unknown file type $type"; return -1}
    }

    return $retval
}

## Reads a json-formatted parameter file
proc CGit::readParamJSON {filename} {

    set fftype ""
    set ffparams ""

    ## Slurp up the file
    set data [read [set fid [open $filename]]][close $fid]

    ## Read the json
    if {[catch {package present "json"}]} {
        package require json
    }
    set data [::json::json2dict $data]

    set fftype [dict get $data "fftype"]
    if { $fftype == {} } {
        cgCon -error "Can't determine force field type"
        cgCon -error "Missing \"fftype\" field?"
        return -code error -1
    }

    set ffdata [dict get $data "params"]
    if { $ffdata == {} } {
        cgCon -error "No parameters defined in file?"
        return -code error -1
    }

    ## Run the sub-parsers
    switch -exact $fftype {
        sdk {set retval [readParamSDK $ffdata]}
        default {cgCon -error "Unknown force field type"}
    }

    return -code ok 0
}

## Reads the SDK-json formatted parameter file
proc CGit::readParamSDK {data} {

    variable params

    set i [lindex [lsort -increasing -integer [dict keys $params]] end]
    if { $i == {} } {set i 0}

    ## Parse the data object by object to give it a unique key
    set nparamsread 0
    foreach obj $data {

        ## Check for exact duplicates
        set id [dict values $params $obj]

        if {$id == ""} {
            dict set params [incr i] $obj
            incr nparamsread
        }
    }

    cgCon -info "Read $nparamsread parameters"
}

# +---------------------------------------+
# |  Analyze the mol specified, assigning |
# |  available topology data to the atoms |
# +---------------------------------------+

## Analyzes the CG molecule, assigns properties, updates bonds, angles, dihedrals, impropers
# \param args molid ?bonds? ?angles? ?dihedrals? ?impropers?
# \return error status
#
# This command applies the properties from the map/topology files to
# the CG particles: mass, charge, index, atom_type segname
#
# If "bonds" and/or "angles" and/or "dihedrals" and/or "impropers" is
# specified, topotools is used to guess those.
#
# \note This command must be run before outputting lammps data/parameter files
proc CGit::reanalyze_mol { args } {

    variable topologies
    variable sys

    set sel [lindex $args 0]
    set seltext [$sel text]
    set flags [lrange $args 1 end]
    set molid [$sel molid]

    ## Default flags
    if {$flags == {}} {set flags {all bonds angles dihedrals}}

    ## Get all atom names from selection
    set unique_atoms [lsort -unique [$sel get {name resname}]]

    ## Assign topology properties to atoms (mass, charge, user, atomtype, segname)
    foreach x $unique_atoms {
        lassign $x name resname

        ## Get properties from the map table
        set resprops [getMap $resname -keys {name mass charge type} -transpose]
        set idx [lsearch -exact -index 0 $resprops $name]

        if { $idx < 0} {
            cgCon -err "missing properties for atom $name for residue $resname"
            return -code error -1
        }

        set mass   [lindex $resprops [list $idx 1]]
        set charge [lindex $resprops [list $idx 2]]
        set type   [lindex $resprops [list $idx 3]]

        ## Get corresponding type-type sigma for radius
        foreach {potential epsilon sigma} {dummy dummy 1.0} break
        set p [getParam pair $type $type -keys {potential epsilon sigma}]
        dict with p {}

        ## Convert to RMIN/2 for LJ 9,6 (rmin = 1.5**(1/3)*sigma)
        set sigma [expr {$sigma * 0.5723}]

        set sel2 [atomselect $molid "(name $name and resname $resname) and ($seltext)"]
        foreach a $flags {
            switch $a {
                none   { break }
                mass   { $sel2 set mass $mass }
                charge { $sel2 set charge $charge }
                type   { $sel2 set type $type }
                radius { $sel2 set radius $sigma}
                all    {
                    $sel2 set mass $mass
                    $sel2 set charge $charge
                    $sel2 set type $type
                    $sel2 set radius $sigma
                }
                default { continue }
            }
        }
        $sel2 delete
    }

    ## Guess bonds, angles, dihedrals, impropers. Skip water.
    set sel2 [atomselect $molid "($seltext) and not water"]

    ## type the bonds, guess the angles etc.., update VMD
    foreach a $flags {
        switch $a {
            bonds {puts "retyping bonds"; topo -sel $sel2 retypebonds}
            angles {puts "retyping angles"; topo -sel $sel2 retypeangles}
            dihedrals {puts "retyping dihedrals"; topo -sel $sel2 retypedihedrals}
            impropers {puts "retyping impropers"; topo -sel $sel2 retypeimpropers}
            default { continue }
        }
    }

    mol reanalyze $molid

    $sel2 delete

    return -code ok $sys(OK)
}


# +------------------------------------------------------------------+
# | Apply topologies read in from the maps to the selected molecules |
# +------------------------------------------------------------------+

proc CGit::set_bonds { args } {

    variable map
    variable sys

    set sel [lindex $args 0]
    set molid [$sel molid]
    set seltext [$sel text]

    ## Delete all the bonds in the selection
    topo -sel $sel clearbonds

    ## Get unique residue types
    set resnames [lsort -unique [$sel get resname]]

    set sdk_bondlist {}
    set topo_bondlist {}
    foreach r $resnames {

        if {[catch {dict get $map $r bonds} sdk_bondlist]} {
            continue
        }

        # Make sure we actually have the bonding topology defined
        # and not just an empty list. If the topology specifies
        # "none"  or "auto" then skip it.
        if { [lsearch $sdk_bondlist "*none*"] >= 0 ||
             [llength $sdk_bondlist] == 0} {continue}

        ## Loop over all atom pairs. This assumes that the selection
        ## above isn't missing any atoms from the residues, and that
        ## the indices are numbered consecutively for each residue.
        
        set pbonds {}
        foreach pair $sdk_bondlist {

            ## Extract residues that make bonds to residues before 
            ## or after (+N / -C) that create a polymer chain  
            if {[lsearch $pair "+*"] >= 0 ||
                [lsearch $pair "-*"] >= 0} {
                  lappend pbonds $pair; continue
            }

            set sel2 [atomselect $molid "(name $pair and resname $r) and ($seltext)"]

            foreach {x y} [$sel2 get index] {
                lappend topo_bondlist [list $x $y]
            }

            $sel2 delete
        }
   
        ## Bonds between residues in the same chain
        foreach pair $pbonds {

          ## Get the names without the +/- identifiers 
          set names [string map  {+ "" - ""} $pair]
         
          ## Residues should be consecutive in chains. Chains are never
          ## connected, this should be done via topotools, if necessary, or
          ## just given the same chain identifier. No checks are currently
          ## done, it is up to the user to make sure the correct topology is
          ## provided and subsequently checked.
          set sel2 [atomselect $molid "(name $names and resname $r) and ($seltext)"] 
          set props [$sel2 get {chain resid index}]
          $sel2 delete

          ## Sort by chain, then by resid
          set props [lsort -index 0 [lsort -index 1 -integer $props]] 

          ## Check for bonds between residues atom indices n-n+1,
          ## n+2-n+3
          foreach {x y} $props {
            lassign $x chainx residuex indexx
            lassign $y chainy residuey indexy
  
            if { $chainx == $chainy &&
                 $residuey == $residuex + 1 } {
                   lappend topo_bondlist [list $indexx $indexy]
            }
          }

          ## Check for bonds between residues atom indices n+1-n+2,
          ## n+3-n+4
          foreach {x y} [lrange $props 1 end] {
            lassign $x chainx residuex indexx
            lassign $y chainy residuey indexy
  
            if { $chainx == $chainy &&
                 $residuey == $residuex + 1 } {
                   lappend topo_bondlist [list $indexx $indexy]
            }
          }
        }
    }

    ## Set bond list and retype
    topo -sel $sel setbondlist $topo_bondlist
    topo -sel $sel retypebonds

    return -code ok $sys(OK)
}

proc CGit::set_angles { args } {

    variable map
    variable sys

    set sel [lindex $args 0]
    set molid [$sel molid]
    set seltext [$sel text]

    ## Delete all the angles in the selection
    topo -sel $sel clearangles

    ## Get unique residue types
    set resnames [lsort -unique [$sel get resname]]

    set sdk_anglelist {}
    set topo_anglelist {}
    foreach r $resnames {

        if {[catch {dict get $map $r angles} sdk_anglelist]} {
            continue
        }

        ## Make sure we actually have the angles defined
        ## and not just an empty list. If the topology
        ## specifies "none" skip it, or guess if "auto"
        if { [lsearch $sdk_anglelist "*none*"] >= 0 ||
             [llength $sdk_anglelist] == 0} {continue}

        ## Guess the angles if we can
        if {  [lsearch $sdk_anglelist "*auto*"] >= 0 } {
            ## Guess angles and store them
            set sel2 [atomselect $molid "($seltext) and resname $r"]
            puts "guessing angles using topotools for residue $r"
            topo -sel $sel2 guessangles
            set topo_anglelist [concat $topo_anglelist\
                                    [topo -sel $sel2 getanglelist]]
            $sel2 delete
        } else {

            ## Loop over all atom tuples.
            set tangles {}
            foreach tuple $sdk_anglelist {

                if {[llength $tuple] != 3} {
                    cgCon -warn "bad angle definition: $tuple"
                }

                ## Extract residues that make angles to residues before 
                ## or after (+N / -C) that create a polymer chain  
                if {[lsearch $tuple "+*"] >= 0 ||
                    [lsearch $tuple "-*"] >= 0} {
                      lappend tangles $tuple; continue
                }

                ## Order here is important! With bonding it doesn't
                ## matter, but we need to be careful and not mess up
                ## the order for angles if the atoms read in by the
                ## user aren't commensurate with the angle list.
                set sel2 [atomselect $molid "(name $tuple and resname $r) and ($seltext)"]
                set name [$sel2 get name]
                set index [$sel2 get index]

                ## Isolate each atom, loop over in parallel, this assumes that all atoms
                ## are ordered the same in each residue.
                set a1 [lsearch -exact -all $name [lindex $tuple 0]]
                set a2 [lsearch -exact -all $name [lindex $tuple 1]]
                set a3 [lsearch -exact -all $name [lindex $tuple 2]]

                if {[llength $a1] != [llength $a2] ||
                    [llength $a1] != [llength $a3]} {
                    cgCon -err "Missing atoms for angle $tuple for residue type $r or bad atom name"
                }

                foreach id1 $a1 id2 $a2 id3 $a3 {
                    lappend topo_anglelist [list "UNK-UNK-UNK"\
                                                [lindex $index $id1]\
                                                [lindex $index $id2]\
                                                [lindex $index $id3]]
                }

                ## Cleanup selection
                $sel2 delete
            }

            foreach tuple $tangles {

              ## Get the names without the +/- identifiers 
              set names [string map  {+ "" - ""} $tuple]
             
              ## Select the atoms
              set sel2 [atomselect $molid "(name $names and resname $r) and ($seltext)"]
              set props [$sel2 get {name chain resid index}]
              $sel2 delete

              ## Sort by chain, then by resid
              set props [lsort -index 1 [lsort -integer -index 2 $props]] 

              ## Isolate each atom, loop over in parallel, this assumes that all atoms
              ## are ordered the same in each residue and that there are no
              ## missing atoms. 
              set a1 [lsearch -inline -exact -all -index 0 $props [lindex $names 0]]
              set a2 [lsearch -inline -exact -all -index 0 $props [lindex $names 1]]
              set a3 [lsearch -inline -exact -all -index 0 $props [lindex $names 2]]

              #^-[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9] -*-*-*
              #^-[A-Za-z0-9]+--[A-Za-z0-9]+-[A-Za-z0-9] -*--*-*
              # ^[A-Za-z0-9]+-\+[A-Za-z0-9]+-\+[A-Za-z0-9] *-+*-+*
              # ^[A-Za-z0-9]+-[A-Za-z0-9]+-\+[A-Za-z0-9] *-*-+*
         
              ## Validate tuple format and adjust lists to match
              ## correct traversal order 
              switch -regexp [join $tuple "-"] {
 
                {^-[A-Za-z0-9]+--[A-Za-z0-9]+-[A-Za-z0-9]+} {
                  set a3 [lrange $a3 1 end]
                }
                
                {^-[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+} {
                  set a2 [lrange $a2 1 end]
                  set a3 [lrange $a3 1 end]
                }
                
                {^[A-Za-z0-9]+-[A-Za-z0-9]+-\+[A-Za-z0-9]+} {
                  set a3 [lrange $a3 1 end]
                }
    
                {^[A-Za-z0-9]+-\+[A-Za-z0-9]+-\+[A-Za-z0-9]+} {
                  set a2 [lrange $a2 1 end]
                  set a3 [lrange $a3 1 end]
                }
                
                default { 
                  cgCon -error "Mangled angle definition: $tuple"
                  return -code error $sys(ERROR)
                }
              } 
          
              foreach id1 $a1 id2 $a2 id3 $a3 {
             
                ## Check for the end of the id list 
                if {$id1 == "" || $id2 == "" || $id3 == ""} {break}

                lassign $id1 name1 chain1 resid1 index1
                lassign $id2 name2 chain2 resid2 index2
                lassign $id3 name3 chain3 resid3 index3

                ## Check for same chain and consecutive resids.
                ## It's messy, but it checks all possible cases
                ## to make sure the topology makes sense
                if {$chain1 == $chain2 && 
                    $chain1 == $chain3 && (
                    ($resid1 == $resid2 && $resid1 + 1 == $resid3) || 
                    ($resid1 + 1 == $resid2 && $resid1 + 1 == $resid3)
                    )} { 
                      lappend topo_anglelist [list "UNK-UNK-UNK"\
                        $index1 $index2 $index3]
                }
              }
           }
        }
     }

    ## Set angle list and retype
    topo -sel $sel setanglelist $topo_anglelist
    topo -sel $sel retypeangles

    return -code ok $sys(OK)
}

proc CGit::set_dihedrals { args } {

    variable map
    variable sys

    set sel [lindex $args 0]
    set molid [$sel molid]
    set seltext [$sel text]

    ## Delete all the dihedrals in the selection
    topo -sel $sel cleardihedrals

    ## Get unique residue types
    set resnames [lsort -unique [$sel get resname]]

    set sdk_dihedrallist {}
    set topo_dihedrallist {}
    foreach r $resnames {

        if {[catch {dict get $map $r dihedrals} sdk_dihedrallist]} {
            continue
        }

        ## Make sure we actually have the dihedrals defined
        ## and not just an empty list. If the topology
        ## specifies "none" skip it, or guess if "auto"
        if { [lsearch $sdk_dihedrallist "*none*"] >= 0 ||
             [llength $sdk_dihedrallist] == 0} {continue}

        ## Guess the dihedrals if we can
        if {  [lsearch $sdk_dihedrallist "*auto*"] >= 0 } {
            ## Guess dihedrals and store them
            set sel2 [atomselect $molid "($seltext) and resname $r"]
            puts "guessing dihedrals using topotools for residue $r"
            topo -sel $sel2 guessdihedrals
            set topo_dihedrallist [concat $topo_dihedrallist\
                                       [topo -sel $sel2 getdihedrallist]]
            $sel2 delete
        } else {

            ## Loop over all atom tuples.
            set tdihed {}
            foreach tuple $sdk_dihedrallist {

                if {[llength $tuple] != 4} {
                    cgCon -warn "bad dihedral definition: $tuple"
                }

                ## Extract residues that make dihedrals to residues before 
                ## or after (+N / -C) that create a polymer chain  
                if {[lsearch $tuple "+*"] >= 0 ||
                    [lsearch $tuple "-*"] >= 0} {
                      lappend tdihed $tuple; continue
                }

                ## Order here is important! With bonding it doesn't
                ## matter, but we need to be careful and not mess up
                ## the order for dihedrals if the atoms read in by the
                ## user aren't commensurate with the dihedral list.
                set sel2 [atomselect $molid "(name $tuple and resname $r) and ($seltext)"]
                set name [$sel2 get name]
                set index [$sel2 get index]

                ## Isolate each atom, loop over in parallel, this assumes that all atoms
                ## are order the same in each residue.
                set a1 [lsearch -exact -all $name [lindex $tuple 0]]
                set a2 [lsearch -exact -all $name [lindex $tuple 1]]
                set a3 [lsearch -exact -all $name [lindex $tuple 2]]
                set a4 [lsearch -exact -all $name [lindex $tuple 3]]

                if {[llength $a1] != [llength $a2] ||
                    [llength $a1] != [llength $a3] ||
                    [llength $a1] != [llength $a4]} {
                    cgCon -err "Missing atoms for dihedral $tuple for residue type $r or bad atom name"
                }

                foreach id1 $a1 id2 $a2 id3 $a3 id4 $a4 {
                    lappend topo_dihedrallist [list "UNK-UNK-UNK-UNK"\
                                                   [lindex $index $id1]\
                                                   [lindex $index $id2]\
                                                   [lindex $index $id3]\
                                                   [lindex $index $id4]]
                }

                ## Cleanup selection
                $sel2 delete
            }
           
            foreach tuple $tdihed {

              ## Get the names without the +/- identifiers 
              set names [string map  {+ "" - ""} $tuple]
             
              ## Select the atoms
              set sel2 [atomselect $molid "(name $names and resname $r) and ($seltext)"]
              set props [$sel2 get {name chain resid index}]
              $sel2 delete

              ## Sort by chain, then by resid
              set props [lsort -index 1 [lsort -integer -index 2 $props]] 

              ## Isolate each atom, loop over in parallel, this assumes that all atoms
              ## are ordered the same in each residue and that there are no
              ## missing atoms. 
              set a1 [lsearch -inline -exact -all -index 0 $props [lindex $names 0]]
              set a2 [lsearch -inline -exact -all -index 0 $props [lindex $names 1]]
              set a3 [lsearch -inline -exact -all -index 0 $props [lindex $names 2]]
              set a4 [lsearch -inline -exact -all -index 0 $props [lindex $names 3]]

              #^[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+-\+[A-Za-z0-9]+ *-*-*-+*
              #^[A-Za-z0-9]+-[A-Za-z0-9]+-\+[A-Za-z0-9]+-\+[A-Za-z0-9]+ *-*-+*-+*
              #^[A-Za-z0-9]+-\+[A-Za-z0-9]+-\+[A-Za-z0-9]+-\+[A-Za-z0-9]+ *-+*-+*-+*
              #^[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+--[A-Za-z0-9]+ *-*-*--*        
              #^[A-Za-z0-9]+-[A-Za-z0-9]+--[A-Za-z0-9]+--[A-Za-z0-9]+ *-*--*--*
              #^[A-Za-z0-9]+--[A-Za-z0-9]+--[A-Za-z0-9]+--[A-Za-z0-9]+ *--*--*--*
 
              ## Validate tuple format and adjust lists to match
              ## correct traversal order 
              switch -regexp [join $tuple "-"] {
 
                {^[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+-\+[A-Za-z0-9]+}  {
                  set a4 [lrange $a4 1 end]
                }
                
                {^[A-Za-z0-9]+-[A-Za-z0-9]+-\+[A-Za-z0-9]+-\+[A-Za-z0-9]+}  {
                  set a3 [lrange $a3 1 end]
                  set a4 [lrange $a4 1 end]
                }
                
                {^[A-Za-z0-9]+-\+[A-Za-z0-9]+-\+[A-Za-z0-9]+-\+[A-Za-z0-9]+}  {
                  set a2 [lrange $a2 1 end]
                  set a3 [lrange $a3 1 end]
                  set a4 [lrange $a4 1 end]
                }
    
                {^[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+--[A-Za-z0-9]+} {
                  set a4 [lrange $a4 1 end]
                }

                {^[A-Za-z0-9]+-[A-Za-z0-9]+--[A-Za-z0-9]+--[A-Za-z0-9]+} {
                  set a3 [lrange $a3 1 end]
                  set a4 [lrange $a4 1 end]
                }

                {^[A-Za-z0-9]+--[A-Za-z0-9]+--[A-Za-z0-9]+--[A-Za-z0-9]+} {
                  set a2 [lrange $a2 1 end]
                  set a3 [lrange $a3 1 end]
                  set a4 [lrange $a4 1 end]
                }
 
                default { 
                  cgCon -error "Mangled dihedral definition: $tuple"
                  return -code error $sys(ERROR)
                }
              } 
          
              foreach id1 $a1 id2 $a2 id3 $a3 id4 $a4 {
             
                ## Check for the end of the id list 
                if {$id1 == "" || $id2 == "" ||
                    $id3 == "" || $id4 == ""} {break}

                lassign $id1 name1 chain1 resid1 index1
                lassign $id2 name2 chain2 resid2 index2
                lassign $id3 name3 chain3 resid3 index3
                lassign $id4 name4 chain4 resid4 index4

                ## Check for same chain and consecutive resids.
                ## It's messy, but it checks all possible cases
                ## to make sure the topology makes sense
                if {$chain1 == $chain2 && 
                    $chain1 == $chain3 &&
                    $chain1 == $chain4 && (
                    ($resid1 + 1 == $resid2 && 
                     $resid1 + 1 == $resid3 && 
                     $resid1 + 1 == $resid4) || 
                    ($resid1 == $resid2 && 
                     $resid1 + 1 == $resid3 && 
                     $resid1 + 1 == $resid4) ||
                    ($resid1 == $resid2 && 
                     $resid1 == $resid3 && 
                     $resid1 + 1 == $resid4)
                    )} { 
                      lappend topo_dihedrallist [list "UNK-UNK-UNK-UNK"\
                        $index1 $index2 $index3 $index4]
                }
              }
           } 
        }
    }

    ## Set dihedral list and retype
    topo -sel $sel setdihedrallist $topo_dihedrallist
    topo -sel $sel retypedihedrals

    return -code ok $sys(OK)
}

## Print out a nice formatted database of values in plain-text
proc CGit::printParam {{fname {stdout}}} {

    variable sys
    variable params

    set fid 0
    if {$fname != "stdout"} {
        if {[catch {open $fname w} fid]} {
            cgCon -err "print_params: Unable to open $fname for writing"
            return -code error $sys(ERROR)
        }

        puts "Writing output to $fname"

    } else {
        set fid $fname
    }

    ## Loop over unique param types
    foreach p [lsort -unique\
                   [getParam * -values -all -keys {param} -uindex 1]] {
        dict for {key value} [getParam $p -all] {
            set value [join [dict values $value]]
            set n [llength $value]
            set fmt [lrepeat $n "%8s"]
            puts $fid [format $fmt {*}$value]
        }
    }

    if {$fname != "stdout"} {close $fid}

    return -code ok
}


# +-----------------+
# | Getters/Setters |
# +-----------------+

## Gets fields based on potential type, returns all values or
## values selected by user. This is the main proc for interacting
## with the potentials.

# set params [getParam pair W W -keys {potential sigma epsilon} -values]
# dict with params {}
# puts $potential ->> returns looked up potential

proc CGit::getParam {args} {

    variable params

    ## Flags By default we return a formatted dict
    set f_all     0
    set f_keys    0; set keys "all"
    set f_indices 0
    set f_values  0
    set f_unique  0
    set f_index   0; set uidx 0

    set nargs [llength $args]
    set newargs ""

    ## Check for flags and remove them from arguments
    for {set i 0} {$i < $nargs} {incr i} {
        set arg [lindex $args $i]
        if {[string match -?* $arg]} {
            set val [lindex $args [expr $i+1]]
            switch -- $arg {
                -all     {set f_all 1}
                -keys    {set f_keys 1; set keys $val; incr i}
                -indices {set f_indices 1}
                -values  {set f_values 1}
                -both    {set f_indices 1; set f_values 1}
                -unique  {set f_unique 1}
                -uindex  {set f_index 1; set uidx $val; incr i}

                default {cgCon -info "Unknown option $arg"; return}
            }
        } else {
            lappend newargs $arg
        }
    }

    set retval ""
    set command [lindex $newargs 0]
    set newargs [lrange $newargs 1 end]

    switch -exact $command {
        bond {

            lassign $newargs t1 t2
            if {$t1 == {}} {set t1 "*"}
            if {$t2 == {}} {set t2 "*"}

            set retval [dict filter $params value\
                            "*bond*\{$t1 $t2\}*"]
            if {$retval == {}} {
                set retval [dict filter $params value\
                                "*bond*\{$t2 $t1\}*"]
            }
        }

        angle {

            lassign $newargs t1 t2 t3
            if {$t1 == {}} {set t1 "*"}
            if {$t2 == {}} {set t2 "*"}
            if {$t3 == {}} {set t3 "*"}

            set retval [dict filter $params value\
                            "*angle*\{$t1 $t2 $t3\}*"]
            if {$retval == {}} {
                set retval [dict filter $params value\
                                "*angle*\{$t3 $t2 $t1\}*"]
            }
        }

        dihed -
        dihedral {

            lassign $newargs t1 t2 t3 t4
            if {$t1 == {}} {set t1 "*"}
            if {$t2 == {}} {set t2 "*"}
            if {$t3 == {}} {set t3 "*"}
            if {$t4 == {}} {set t4 "*"}

            set retval [dict filter $params value\
                            "*dihedral*\{$t1 $t2 $t3 $t4\}*"]
            if {$retval == {}} {
                set retval [dict filter $params value\
                                "*dihedral*\{$t4 $t3 $t2 $t1\}*"]
            }
        }

        lj -
        pair {

            lassign $newargs t1 t2
            if {$t1 == {}} {set t1 "*"}
            if {$t2 == {}} {set t2 "*"}

            set retval [dict filter $params value\
                            "*pair*\{$t1 $t2\}*"]
            if {$retval == {}} {
                set retval [dict filter $params value\
                                "*pair*\{$t2 $t1\}*"]
            }
        }

        default {## This can be abused...
            set retval [dict filter $params value $command] }
    }

    ## Return a catenated list of the {key value} pairs
    ## requested
    if {$f_keys} {
        set newretval {}
        dict for {key val} $retval {
            set innerlist {}
            foreach k $keys {
                lappend innerlist [dict filter $val key $k]
            }
            lappend newretval $key [join $innerlist]
        }
        set retval $newretval
        unset newretval
    }

    ## Return just keys or values?
    if {$f_indices && !$f_values} {
        return [lsort -integer -increasing\
                    [dict keys $retval *]]
    } elseif {$f_values && !$f_indices} {
        set retval [dict values $retval]
    }

    ## Get unique values based on specified index
    ## We need all 3 flags for this to make sense
    if {$f_unique && $f_index && $f_values} {
        #set retval [lsearch -all -inline -index $uidx -subindices $retval *]
        #set retval [lsort -unique $retval]
        set retval [lsort -unique -index $uidx $retval]
    } elseif {$f_index && $f_values} {
        set retval [lsearch -all -inline -index $uidx -subindices $retval *]
    }

    ## Return all matching values if requested
    if {$f_all} {return $retval}

    ## Return the dict for the matching value
    return [lindex $retval end]
}

## Generate pair parameters via selected combination rules
proc CGit::pairCombine {atype1 atype2 {rule geometric}} {

    set retval {}

    ## Check if there exists a usable combination already definied
    foreach {potential epsilon sigma} {"" "" ""} break
    set p [getParam pair $atype1 $atype2 -keys {potential epsilon sigma}]
    if {$p != ""} {return $p}

    foreach {potential epsilon sigma} {"" "" ""} break
    set p [getParam pair $atype1 $atype1 -keys {potential epsilon sigma}]
    dict with p {
        set pottype1 $potential; set epsilon1 $epsilon; set sigma1 $sigma
    }

    foreach {potential epsilon sigma} {"" "" ""} break
    set p [getParam pair $atype2 $atype2 -keys {potential epsilon sigma}]
    dict with p {
        set pottype2 $potential; set epsilon2 $epsilon; set sigma2 $sigma
    }

    if {$epsilon1 == "" || $epsilon2 == "" ||
        $sigma1 == "" || $sigma2 == ""} { return {} }

    ## Check if we have a water interaction
    if {[string match $pottype1 "lj12_4"] ||
        [string match $pottype2 "lj12_4"]} {
        set pottype "lj12_4"
    } else {
        set pottype "lj9_6"
    }

    ## Return the correct interaction type and calculated values
    switch -exact $rule {
        geometric {
            set s [expr {sqrt(double($sigma1) * double($sigma2))}]
            set e [expr {sqrt(double($epsilon1) * double($epsilon2))}]
        }
        arithmetic {
            set s [expr {0.5 * (double($sigma1) + double($sigma2))}]
            set e [expr {0.5 * (double($epsilon1) + double($epsilon2))}]
        }
        lb {
            set s [expr {0.5 * (double($sigma1) + double($sigma2))}]
            set e [expr {sqrt(double($epsilon1) * double($epsilon2))}]
        }
        default {
            set s 0.0; set e 0.0;
        }
    }

    return [list potential $pottype epsilon $e sigma $s]
}

## Read param wrapper based on provided or guessed filetype
proc CGit::readTopo {filename {type guess}} {

    cgCon -info "Reading topologies from $filename"

    if {![string compare $type "guess"]} {
        set type [file extension $filename]
    }

    switch -exact $type {
        ".json" -
        "json"  {set retval [readTopoJSON $filename]}
        default {cgCon -err "Unknown file type $type"; return -1}
    }

    return $retval
}

## Reads a json-formatted parameter file
proc CGit::readTopoJSON {filename} {

    variable map

    set fftype ""
    set ffparams ""

    ## Slurp up the file
    set data [read [set fid [open $filename]]][close $fid]

    ## Read the json
    if {[catch {package present "json"}]} {
        package require json
    }
    set data [::json::json2dict $data]

    set fftype [dict get $data "fftype"]
    if { $fftype == {} } {
        cgCon -error "Can't determine force field type"
        cgCon -error "Missing \"fftype\" field?"
        return -code error -1
    }

    set ffdata [dict get $data "topo"]
    if { $ffdata == {} } {
        cgCon -error "Missing 'topo' keyword, no topologies defined in file?"
        return -code error -1
    }

    ## Merge in the dictionary
    switch -exact $fftype {
        sdk {set map [dict merge $map $ffdata]}
        default {cgCon -error "Unknown force field type"}
    }

    cgCon -info "Topology database contains [dict size $map] entries"

    return -code ok 0
}

## Utility to get map properties at the coast of some
## overhead for parsing the flags and transposing
## if necessary
proc ::CGit::getMap {args} {

    variable map

    ## Flags
    set f_keys    0; set keys {type name map charge mass}
    set f_trans   0;

    set nargs [llength $args]
    set newargs ""

    ## Check for flags and remove them from arguments
    for {set i 0} {$i < $nargs} {incr i} {
        set arg [lindex $args $i]
        if {[string match -?* $arg]} {
            set val [lindex $args [expr $i+1]]
            switch -- $arg {
                -keys    {set f_keys 1; set keys $val; incr i}
                -transpose {set f_trans 1}
                default {cgCon -info "Unknown option $arg"; return}
            }
        } else {
            lappend newargs $arg
        }
    }

    set retval ""
    set resname [lindex $newargs 0]

    ## Speed things up a bit if we're requesting 1 thing
    if {[llength $keys] == 1} {

        if {[catch {dict get $map $resname $keys} retval]} {
            cgCon -err "Unknown property \"$keys\" for residue \"$resname\""
            return -code error ""
        }

    } else {

        ## Get each property and construct sublists
        set retval {}
        foreach k $keys {
            if {[catch {dict get $map $resname $k} val]} {
                cgCon -err "Unknown property \"$k\" for residue \"$resname\""
                return -code error ""
            }
            lappend retval $val
        }

        ## Transpose the properties to return "by bead"
        ## so that the string can be applied directly
        ## to an atomselection

        if {$f_trans} {set retval [transpose $retval]}
    }

    return $retval
}


## lsearch hack to transpose matrix
proc ::CGit::transpose {matrix} {
    set res {}
    for {set index 0} {$index < [llength [lindex $matrix 0]]} {incr index} {
        if {[catch {lsearch -all -inline -subindices\
                        -index $index $matrix *} val]} {
            cgCon -error "Unequal number of properties for bead"
            return -code error
        }
        lappend res $val
    }
    return $res
}

## For the relatively small arrays we're transposing
## this method wins out by a few microseconds per iteration
## than the lsearch transpose hack above, but this lacks implicit
## error handling...
proc ::CGit::transpose2 {matrix} {
    set res {}
    for {set j 0} {$j < [llength [lindex $matrix 0]]} {incr j} {
        set newrow {}
        foreach oldrow $matrix {
            lappend newrow [lindex $oldrow $j]
        }
        lappend res $newrow
    }
    return $res
}


# +--------------------+
# | Cleanup Everything |
# +--------------------+

proc CGit::reset {} {
    ## Clear everything
    clear all
}

# +--------------------+
# |  Clear Everything  |
# +--------------------+

proc CGit::clear {{flag all}} {

    variable map
    variable params

    switch -- $flag {

        map -
        maps {
            cgCon -info "Clearing maps/topologies"
            set map ""
        }

        param -
        params {
            cgCon -info "Clearing parameters"
            set params {}
        }

        all {
            ## Clear parameters
            cgCon -info "Clearing parameters"
            set params {}

            ## Clear the maps/topologies
            cgCon -info "Clearing maps/topologies"
            set map ""
        }
        default {}
    }
}
