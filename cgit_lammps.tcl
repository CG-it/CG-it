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

# +----------------------------------------------------------+
# | just a wraper for topo writelammpsdata and animate write |
# +----------------------------------------------------------+

proc CGit::writelammpsdata {molid fname {style full}} {
    topo -molid $molid writelammpsdata $fname $style
}

proc CGit::readlammpsdata {molid fname {style full}} {

  ## Read the topology file and do the best we can assigning 
  ## atom names to the molecules. There are a lot of assumptions:
  # 0. There are resnames to begin with
  # 1. Residues are grouped contigiously w.r.t the atom selection 
  # 2. The atom order in each residue type is consistent across residues 
  # 3. There are no missing atoms
  # 4. There is not an insane (10e6) amount of beads/atoms 
  
  set newmolid [topo -molid $molid readlammpsdata $fname $style]

  ## Get list of residues, and their numbers
  set sel      [atomselect $newmolid "all"] 
  set residues [lsort -unique [$sel get resname] 
  $sel delete

  array set nres {}
  array set sels {}
  
  ## loop over residue types, calculate the number of residues each
  foreach r $residues {
    set sel [atomselect $newmolid "resname $r"]
    lappend sels($r) $sel
    lappend nres($r) nres $sel
  }

  ## create a list with atom names * nres
  set names {} 
  foreach r $residues {
    set names [::CGit::getMap $r -keys name]
    set names [lrepeat $nres($r) $names] 

    ## Things get dicey here. If there are too many atoms
    ## to apply here via the atomselections, vmd could puke
    ## if we overflow the selection text buffer 
    $sels($r) set name $names
  }

  #cleanup 
  foreach r $residues { $sels($r) delete }

  return $newmolid
}

# +---------------------------------------------+
# | Write out a lammps formatted parameter file |
# +---------------------------------------------+

proc CGit::writelammpsparam {molid fname {flag none} {guess 0}} {

    variable params

    ## Conveniently, Topotools enumerates
    ## all the bond, angle, dihedral types etc..
    ## We just need to assign indices to the types
    ## and then enumerate over all possible combinations
    ## for the non-bonded interactions and write out the others
    ## directly into a format compatible with LAMMPS.

    if {[catch {open $fname "w"} fid]} {
        cgCon -err "Can't open file $fname for writing"
        return
    }

    ## Write a preamble
    lammps_preamble $fid $flag

    ## Get atom properties
    set sel [atomselect $molid "all"]
    set atom_props [lsort -unique -ascii -index 0 [$sel get {type mass charge}]]

    set types      [lsearch -all -index 0 -subindices -inline $atom_props *]
    set masses     [lsearch -all -index 1 -subindices -inline $atom_props *]
    set charges    [lsearch -all -index 2 -subindices -inline $atom_props *]

    $sel delete

    ## Convert Bond/angle/dihedral lists to an indexable dictionary
    dict set topodict bond [topo -molid $molid bondtypenames]
    dict set topodict angle [topo -molid $molid angletypenames]
    dict set topodict dihedral [topo -molid $molid dihedraltypenames]

    ## Lookup table to convert between types and their associated numerical id
    ## and back again. Types are number 1 to N..
    catch {array unset lookup}
    foreach x $types {
        set lookup($x) [lsearch -sorted -ascii $types $x]
        incr lookup($x)

        ## Reverse lookup entry
        set lookup($lookup($x)) $x
    }

    puts $fid ""

    ## Masses
    foreach x $types y $masses {
        puts $fid [format "mass  %4.0f  %8.4f \# %s" $lookup($x) $y $x]
    }

    puts $fid ""

    #Non-Bonding Interations

    ## Charged atoms for pair hybrid/overlay with coul/long
    set charged_ids [lsearch -all -not -real $charges 0.0]
    set ci {}
    foreach x $charged_ids {
        lappend ci [expr {$x + 1}]
    }
    set charged_ids $ci
    unset ci

    set i 0
    foreach id1 $charged_ids {
        foreach id2 [lrange $charged_ids $i end] {
            puts $fid [format "\#pair_coeff  %4.0f  %4.0f coul/long \# %s %s"\
                           $id1 $id2 $lookup($id1) $lookup($id2)]
        }
        incr i
    }

    ## Required for overlay and neighborlist multi
    if {[llength $charged_ids] > 0} {
        puts $fid "\#group charged type $charged_ids"
        puts $fid "\#atom_modify first charged"
    }

    puts $fid ""

    ## LJ/SDK for all
    set i 0
    foreach atype1 $types {
        foreach atype2 [lrange $types $i end] {

            foreach {potential epsilon sigma} {"" "" ""} break
            set p [getParam pair $atype1 $atype2 -keys {potential epsilon sigma}]
            dict with p {}

            ## Check if we have any missing parameters
            if {$potential == "" || $epsilon == "" || $sigma == ""} {
		if {$guess} {
		    cgCon -warn "Assigning parameters for $atype1 $atype2 via combining rules"
		    set p [pairCombine $atype1 $atype2 geometric]; dict with p {}
		} else {
		    cgCon -warn "Missing non-bonding parameters for $atype1 $atype2"
		    continue
		}
            }

            puts $fid [format "pair_coeff  %4.0f  %4.0f  %7s  %8.4f  %8.4f \# %s %s"\
                           $lookup($atype1) $lookup($atype2) $potential\
                           $epsilon $sigma $atype1 $atype2]
        }
        incr i
    }

    puts $fid "\n"
    flush $fid

    #if {0} {writecsv $types}

    ## bonded interactions
    set bondtypelist [dict get $topodict bond]
    if {$bondtypelist == {} } {
        vmdcon -warn "No bondtypes found."
    } else {

        foreach x $bondtypelist {
            lassign [split $x "-"] atype1 atype2

            foreach {potential k r0} {"" "" ""} break
            set p [getParam bond $atype1 $atype2 -keys {potential k r0}]
            dict with p {}

            ## Check if we have any missing parameters
            if {$potential == "" || $k == "" || $r0 == ""} {
                cgCon -warn "Missing bonding parameters for $atype1 $atype2"
                continue
            }

            set type [lsearch -ascii $bondtypelist $x]
            incr type;

            puts $fid [format "bond_coeff  %4.0f  %8.4f  %8.4f \# %s-%s"\
                           $type $k $r0 $atype1 $atype2]

        }

        puts $fid "\n"
        flush $fid
    }

    ## Angle types
    set angletypelist [dict get $topodict angle]

    if {$angletypelist == {} } {
        vmdcon -warn "No angletypes found."
    } else {

        ## Check if we have a unique number of angle types, if now, do a hybrid
        set nangtypes [llength [lsort -unique [getParam angle -values -all -keys {potential} -uindex 1]]]
        cgCon -info "System has $nangtypes unique angle potentials"

        foreach x $angletypelist {

            lassign [split $x "-"] atype1 atype2 atype3

            foreach {potential k theta0} {"" "" ""} break
            set p [getParam angle $atype1 $atype2 $atype3 -keys {potential k theta0}]
            dict with p {}
            set angpotential $potential

            ## Check if we have any missing parameters
            if {$potential == "" || $k == "" || $theta0 == ""} {
                cgCon -warn "Missing angle parameters for $atype1 $atype2 $atype3"
                continue
            }

            ## Non-bonding interactions between atoms 1 and 3
            foreach {potential epsilon sigma} {"" "" ""} break
            set p [getParam pair $atype1 $atype3 -keys {potential epsilon sigma}]
            dict with p {}

            ## Check if we have any missing parameters
            if {$potential == "" || $epsilon == "" || $sigma == ""} {
		if {$guess} {
		    cgCon -warn "Assigning parameters for $atype1 $atype3 via combining rules"
		    set p [pairCombine $atype1 $atype3 geometric]; dict with p {}
		} else {
		    cgCon -warn "Missing non-bonding parameters for $atype1 $atype2"
		    continue
		}
            }

            set type [lsearch -ascii $angletypelist $x]
            incr type

            if {$nangtypes > 1} {

                switch -exact $angpotential {
                    harm -
                    harmonic { # Harmonic + Hybrid
                        puts $fid [format "angle_coeff  %4.0f %8s  %8.4f  %8.4f \# %s-%s-%s"\
                                       $type $angpotential $k $theta0\
                                       $atype1 $atype2 $atype3]
                    }

                    sdk { # SDK + Hybrid
                        puts $fid [format "angle_coeff  %4.0f %8s  %8.4f  %8.4f  %7s  %8.4f  %8.4f \# %s-%s-%s"\
                                       $type $angpotential $k $theta0 $potential\
                                       $epsilon $sigma $atype1 $atype2 $atype3]

                    }
                    default {}
                }

            } else { ## Assume non-hybrid sdk angle potential
                puts $fid [format "angle_coeff  %4.0f  %8.4f  %8.4f  %7s  %8.4f  %8.4f \# %s-%s-%s"\
                               $type $k $theta0 $potential $epsilon\
                               $sigma $atype1 $atype2 $atype3]
            }
        }

        puts $fid "\n"
        flush $fid
    }

    ## Dihedral types
    set dihedraltypelist [dict get $topodict dihedral]

    if {$dihedraltypelist == {} } {
        vmdcon -warn "No dihedraltypes found."
    } else {

        foreach x $dihedraltypelist {

            lassign [split $x "-"] atype1 atype2 atype3 atype4

            foreach {potential k n d} {"" "" "" "" ""} break
            set p [getParam dihedral $atype1 $atype2 $atype3 -keys {potential k n d}]
            dict with p {}

            ## Check if we have any missing parameters
            if {$potential == "" || $k == "" || $n == "" || $d == ""} {
                cgCon -warn "Missing dihedral parameters for $atype1 $atype2 $atype3 $atype4"
                continue
            }

            set type [lsearch -ascii $dihedraltypelist $x]
            incr type;

            puts $fid [format "dihedral_coeff  %4.0f  %8.4f  %2d  %5d  0.0 \# %s-%s-%s-%s"\
                           $type $k $n $d \
                           $atype1 $atype2 $atype3 $atype4]
        }
    }

    ## Impropers not supported yet

    close $fid

    return -code ok
}

proc CGit::lammps_preamble {fid {flag none}} {

    ## Useful style configurations

    ## Things that will always be printed
    puts $fid "\# Generated by CGTools - writelammpsparam on [clock format [clock seconds]]"

    ## Optional things
    switch $flag {

        cpu {
            puts $fid "\#pair_style      lj/sdk                    15.0"
            puts $fid "\#pair_style      lj/sdk/coul/long          15.0"
            puts $fid "\#pair_style      hybrid/overlay  lj/sdk 15.0  coul/long 25.0"
        }

        gpu {
            puts $fid "\#pair_style      lj/sdk/coul/long/gpu      15.0"
            puts $fid "\#pair_style      lj/sdk/gpu                15.0"
        }

        both {
            puts $fid "\#pair_style      lj/sdk/coul/long          15.0"
            puts $fid "\#pair_style      lj/sdk                    15.0"
            puts $fid "\#pair_style      hybrid/overlay  lj/sdk 15.0  coul/long 25.0"
            puts $fid "\#pair_style      lj/sdk/gpu                15.0"
            puts $fid "\#pair_style      lj/sdk/coul/long/gpu      15.0"
        }

        none { return }

        default { return }
    }

    ## Model specific settings

    puts $fid "\#bond_style      harmonic"
    puts $fid "\#angle_style     sdk"
    puts $fid "\#angle_style     hybrid sdk harmonic"
    puts $fid "\#dihedral_style  charmm"
    puts $fid "\#special_bonds   lj/coul 0.0 0.0 1.0\n"
}


## Write out a csv for the vdW interaction table
proc CGit::writecsv {types} {

    set fid [open types.csv "w"]

    ## LJ/SDK for all
    set i 0
    foreach atype1 $types {
        foreach atype2 [lrange $types $i end] {
            foreach {potential epsilon sigma} {"" "" ""} break
            set p [getParam pair $atype1 $atype2 -keys {potential epsilon sigma}]
            dict with p {}	    
	    puts $fid [format "%s,%s,%s,%.4f,%.4f"\
			   $atype1 $atype2 $potential $epsilon $sigma]
	}
	incr i
    }
    close $fid
}
