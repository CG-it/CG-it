## Legacy Routines for reading/writing/getting topology and parameter
## data

namespace eval ::CGit::Legacy {
    namespace export Legacy
}

## Reads the Legacy SDK topology file
# \param infile Topology file to be read
# \return error status
#
# This command reads the slightly-modified legacy SDK topology files
# having the following format:
# RESI residue_name
#
# atom index residue_name atom_name atom_type mass charge comment
# atom index residue_name atom_name atom_type mass charge comment
# atom ....
#
# bond 1 2
# bond 1 3
# bond ....
#
# angle 1 2 3
# angle 2 3 4
# angle ...
#
# dihedral 1 2 3 4
# dihedral 2 3 4 5
# dihedral ...
#
#\note Topology files have been depreciated in favor of "map" files.
#\note but we continue their support for the time being
proc CGit::Legacy::readTop {infile {arr null}} {

    variable sys
    upvar $arr topo

    if {[catch {open $infile "r"} fid]} {
        cgCon -err "Can't open topology file $infile for reading"
        return -code error $sys(ERROR)
    }

    cgCon -info "Reading topo from $infile"

    ## Read the topology file line by line
    set resi_read {}
    set label UNK
    set atomi 0
    while {[gets $fid line] >= 0} {

        ## Ignore blanks and comment lines
        if {[llength $line] < 1 || [string match \#* $line]} {continue}

        ## Strip all comments from input
        set line [regsub {\#.*} $line ""]

        ## Make an array "topo" to store the data
        lassign $line desc
        switch -exact -- $desc {

            "RESI" {set label [lindex $line 1]; lappend resi_read $label}
            "atom" {
                lassign $line desc idx resname atomname type mass charge
                lappend topo(idx) $idx
                lappend topo($label,$idx) $atomi
                lappend topo(resname) $label
                lappend topo(name) $atomname
                lappend topo(type) $type
                lappend topo(mass) $mass
                lappend topo(charge) $charge
                incr atomi
            }
            "bond" {
                lassign $line desc idx1 idx2
                lappend topo([list bond $label])\
                    [list $idx1 $idx2]
            }
            "angle" {
                lassign $line desc idx1 idx2 idx3
                lappend topo([list angle $label])\
                    [list $idx1 $idx2 $idx3]
            }
            "dihedral" {
                lassign $line desc idx1 idx2 idx3 idx4
                lappend topo([list dihed $label])\
                    [list $idx1 $idx2 $idx3 $idx4]
            }
            "improper" {
                lassign $line desc idx1 idx2 idx3 idx4
                lappend topo([list impr $label])\
                    [list $idx1 $idx2 $idx3 $idx4]
            }
            "pair" {
                lassign $line desc idx1 idx2
                lappend topo([list pair $label])\
                    [list $idx1 $idx2]
            }
            default {cgCon -warn "Unrecognized topology descriptor $desc"}
        }
    }

    close $fid

    set n [llength $resi_read]
    puts "Read $n residues from file $infile: [join $resi_read]"
    return -code ok $sys(OK)
}

# Reads the Legacy SDK parameter file
# \param infile parameter file to be read
# \return error status
#
# This command reads the legacy SDK parameter file having the
# following format:
#
# bond  atom_type_1 atom_type_2 K (kcal/angstrom**2) r_0 (angstroms)
# angle atom_type_1 atom_type_2 atom_type_3 K (kcal/radian**2) theta_0 (degrees)
# pair  atom_type_1 atom_type_2 lj9_6/lj9_12 epsilon (kcal/mol) sigma (angstrom)
#
proc CGit::Legacy::readParam { args } {

    variable params
    variable sys

    lassign $args infile

    ## Do some globbing if requested
    if {![file isfile $infile]} {
        foreach f [glob $infile] {readParam $f}
        return
    }

    if {[catch {open $infile "r"} fid]} {
        cgCon -err "Can't open paramter file $infile for reading"
        return -code error $sys(ERROR)
    }

    cgCon -info "Reading parameters from $infile"

    ## Read the parameter file line by line
    while {[gets $fid line] >= 0} {

        ## Ignore blanks and comment lines
        if {[llength $line] < 1 || [string match \#* $line]} {continue}

        lassign $line desc

        ## Different key assignments if necessary for particular fields
        switch -exact -- $desc {
            "pair" {

                lassign $line desc atype1 atype2 pottype epsilon sigma

                set key   [list $desc $atype1 $atype2]
                set key2  [list $desc $atype2 $atype1]
                set value [list $pottype $epsilon $sigma]

                ## Check if this parameter is already defined (in reverse order)
                if {[info exists params($key2)]} {
                    set key $key2
                }

                ## Check if this parameter is already defined, and warn the user
                if {[info exists params($key)]} {
                    cgCon -warn "Parameters exist for $key.\n Replacing $params($key) with $value from $infile"
                }

                set params($key) $value
            }

            "bond" {

                lassign $line desc atype1 atype2 k req

                set key   [list $desc $atype1 $atype2]
                set key2  [list $desc $atype2 $atype1] ;# Reverse order
                set value [list $k $req]

                ## Check if this parameter is already defined (in reverse order)
                if {[info exists params($key2)]} {
                    set key $key2
                }

                ## Check if this parameter is already defined, and warn the user
                if {[info exists params($key)]} {
                    cgCon -warn "Parameters exist for $key.\n Replacing $params($key) with $value from $infile"
                }

                set params($key) $value

            }

            "angle" {

                lassign $line desc atype1 atype2 atype3 k thetaeq

                set key   [list $desc $atype1 $atype2 $atype3]
                set key2  [list $desc $atype3 $atype2 $atype1] ;# Reverse Order
                set value [list $k $thetaeq]

                ## Check if this parameter is already defined (in reverse order)
                if {[info exists params($key2)]} {
                    set key $key2
                }

                ## Check if this parameter is already defined, and warn the user
                if {[info exists params($key)]} {
                    cgCon -warn "Parameters exist for $key.\n Replacing $params($key) with $value from $infile"
                }

                set params($key) $value
            }

            "dihedral" {cgCon -warn "dihedrals currently unsupported"}
            "improper" {cgCon -warn "impropers currently unsupported"}

            default {cgCon -warn "Unknown parameter descriptor $desc"}
        }
    }
    return -code ok $sys(OK)
}

## Based on the current molecular connectivity,
## write out bond, angle, dihedral fields in the
## legacy sdk format.
proc CGit::Legacy::writeSDKtopo {sel {fname stdout} {renumber 1}} {

    variable sys

    set molid [$sel molid]
    set seltext [$sel text]

    set fid 0
    if {$fname != "stdout"} {
	if {[catch {open $fname w} fid]} {
	    cgCon -err "writesdktopo: Unable to open $fname for writing"
	    return -code error $sys(ERROR)
	}

	puts "Writing output to $fname"

    } else {
	set fid stdout
    }

    ## Write out records for each residue type in selection
    set resname [lsort -unique -index 0 [$sel get {resname residue}]]
    foreach r $resname {
	lassign $r r rid
	set sel2 [atomselect $molid "($seltext) and residue $rid"]
	if {[$sel2 num] == 0} {
	    continue
	}

	puts $fid "RESI $r"
	puts $fid "\n"

	## Make lookup table for re-numbering indices if requested
	set i 0
	if {$renumber} {
	    foreach index [$sel2 get index] {
		set idxlookup($index) [incr i]
	    }
	} else {
	    foreach index [$sel2 get index] {
		set idxlookup($index) $index
	    }
	}

	## Atom Record
	foreach name [$sel2 get name]\
	    index [$sel2 get index]\
	    resname [$sel2 get resname]\
	    name [$sel2 get name]\
	    type [$sel2 get type]\
	    mass [$sel2 get mass]\
	    charge [$sel2 get charge]\
	    {

		puts $fid [format "atom  %5d%6s%6s%6s%8.3f%8.3f \#%s"\
			       $idxlookup($index) $resname $name $type $mass $charge $resname]
	    }

	puts $fid "\n"

	set atom_props [lsort -unique -ascii -index 0 [$sel2 get {type index}]]
	set atom_props [lsort -increasing -integer -index 1 $atom_props]
	set types      [lsearch -all -index 0 -subindices -inline $atom_props *]
	set index      [lsearch -all -index 1 -subindices -inline $atom_props *]

	## Lookup table to convert between types and their associated numerical id
	catch {array unset lookup}
	foreach x $types idx $index {
	    set lookup($x) $idxlookup($idx)
	}

	#Non-Bonding
	set i 0
	foreach atype1 $types {
	    foreach atype2 [lrange $types $i end] {

		puts $fid [format "pair %5d %5d \#%s-%s" \
			       $lookup($atype1) $lookup($atype2) $atype1 $atype2]

	    }
	    incr i
	}

	puts $fid "\n"

	## Bond Record
	foreach x [lsort -increasing -integer -index 0 [topo -sel $sel2 getbondlist type]] {
	    lassign $x idx1 idx2 type
	    puts $fid [format "bond  %5d %5d \#%s"\
			   $idxlookup($idx1) $idxlookup($idx2) $type]
	}

	puts $fid "\n"

	## Angle Record
	foreach x [lsort -increasing -integer -index 1 [topo -sel $sel2 getanglelist]] {
	    lassign $x name idx1 idx2 idx3
	    puts $fid [format "angle %5d %5d %5d \#%s"\
			   $idxlookup($idx1) $idxlookup($idx2) $idxlookup($idx3) $name]
	}

	puts $fid "\n"

	## Dihedral Record
	foreach x [lsort -increasing -integer -index 1 [topo -sel $sel2 getdihedrallist]] {
	    lassign $x name idx1 idx2 idx3 idx4
	    puts $fid [format "dihedral %5d %5d %5d %5d \#%s"\
			   $idxlookup($idx1) $idxlookup($idx2) $idxlookup($idx3) $idxlookup($idx4) $name]
	}

	$sel2 delete
    }

    if {$fname != "stdout"} {close $fid}

    return -code ok $sys(OK)
}

proc CGit::Legacy::setParam args {

    variable sys
    variable params

    set narg [llength $args]
    lassign $args command

    switch $command {

	bond {

	    lassign $args command t1 t2 k req

	    if {$narg < 3} {
		cgCon -err "setbondparam missing argument: setbondparam <type1> <type2> <k> <req>"
		return -code error $sys(ERROR)
	    }

	    set key  [list bond $t1 $t2]
	    set key2 [list bond $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    ## Check for empty values, use previously defined ones
	    if {[info exists params($key)]} {
		lassign $params($key) k_old req_old
		if {$k == ""} {set k $k_old}
		if {$req == ""} {set req $req_old}

	    }

	    set params($key) [list $k $req]
	}

	angle {

	    lassign $args command t1 t2 t3 k thetaeq

	    if {$narg < 4} {
		cgCon -err "setangleparam missing argument: setangleparam <type1> <type2> <type3> <k> <thetaeq>"
		return -code error $sys(ERROR)
	    }

	    set key  [list angle $t1 $t2 $t3]
	    set key2 [list angle $t3 $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    ## Check for empty values, use previously defined ones if they exist
	    if {[info exists params($key)]} {
		lassign $params($key) k_old thetaeq_old
		if {$k == ""} {set k $k_old}
		if {$thetaeq == ""} {set thetaeq $thetaeq_old}
	    }

	    set params($key) [list $k $thetaeq]
	}

	dihedral {

	    lassign $args command t1 t2 t3 t4 k thetaeq

	    if {$narg < 5} {
		cgCon -err "setdihedralparam missing argument: setdihedralparam <type1> <type2> <type3> <type4> <k> <thetaeq>"
		return -code error $sys(ERROR)
	    }

	    set key  [list dihedral $t1 $t2 $t3 $t4]
	    set key2 [list dihedral $t4 $t3 $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    ## Check for empty values, use previously defined ones if they exist
	    if {[info exists params($key)]} {
		lassign $params($key) k_old thetaeq_old
		if {$k == ""} {set k $k_old}
		if {$thetaeq == ""} {set thetaeq $thetaeq_old}
	    }

	    set params($key) [list $k $thetaeq]
	}

	pair {

	    lassign $args command t1 t2 pottype epsilon sigma

	    if {$narg < 3} {
		cgCon -err "setpairparam missing argument: setpairparam <type1> <type2> <pottype> <epsilon> <sigma>"
		return -code error $sys(ERROR)
	    }

	    set key  [list pair $t1 $t2]
	    set key2 [list pair $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    ## Check for empty values, use previously defined ones if they exist
	    if {[info exists params($key)]} {
		## Check for empty values, use previously defined ones
		lassign $params($key) pottype_old epsilon_old sigma_old
		if {$pottype == ""} {set pottype $pottype_old}
		if {$epsilon == ""} {set epsilon $epsilon_old}
		if {$sigma == ""} {set sigma $sigma_old}
	    }

	    set params($key) [list $pottype $epsilon $sigma]
	}

	default {
	    cgCon -error "set*param, unknown command: $command. bond, angle, dihedral, pair"
	    return -code error $sys(ERROR)
	}
    }

    return -code ok
}

proc CGit::Legacy::deleteParam rgs {

    variable sys
    variable params

    set narg [llength $args]
    lassign $args command

    switch $command {

	bond {

	    lassign $args command t1 t2 k req

	    if {$narg < 3} {
		cgCon -err "deletebondparam missing argument: deletebondparam <type1> <type2> <k> <req>"
		return -code error $sys(ERROR)
	    }

	    set key  [list bond $t1 $t2]
	    set key2 [list bond $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    ## Check for existance in array
	    if {[info exists params($key)]} {
		unset params($key)
	    }
	}

	angle {

	    lassign $args command t1 t2 t3 k thetaeq

	    if {$narg < 4} {
		cgCon -err "deleteangleparam missing argument: deleteangleparam <type1> <type2> <type3> <k> <thetaeq>"
		return -code error $sys(ERROR)
	    }

	    set key  [list angle $t1 $t2 $t3]
	    set key2 [list angle $t3 $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    ## Check for existance in array
	    if {[info exists params($key)]} {
		unset params($key)
	    }
	}

	dihedral {

	    lassign $args command t1 t2 t3 t4 k thetaeq

	    if {$narg < 5} {
		cgCon -err "deletedihedralparam missing argument: deletedihedralparam <type1> <type2> <type3> <type4> <k> <thetaeq>"
		return -code error $sys(ERROR)
	    }

	    set key  [list dihedral $t1 $t2 $t3 $t4]
	    set key2 [list dihedral $t4 $t3 $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    ## Check for existance in array
	    if {[info exists params($key)]} {
		unset params($key)
	    }
	}

	pair {

	    lassign $args command t1 t2 pottype epsilon sigma

	    if {$narg < 3} {
		cgCon -err "deletepairparam missing argument: deletepairparam <type1> <type2> <pottype> <epsilon> <sigma>"
		return -code error $sys(ERROR)
	    }

	    set key  [list pair $t1 $t2]
	    set key2 [list pair $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    ## Check for existance in array
	    if {[info exists params($key)]} {
		unset params($key)
	    }
	}

	default {
	    cgCon -err "delete*param, unknown command: $command. bond, angle, dihedral, pair"
	    return -code error $sys(ERROR)
	}
    }

    return -code ok

}

proc CGit::Legacy::getParam args {

    variable sys
    variable params

    set narg [llength $args]
    lassign $args command

    switch $command {

	bond {

	    if {$narg < 3} {
		cgCon -err "getbondparam missing argument: getbondparam <type1> <type2>"
		return -code ok
	    }

	    lassign $args command t1 t2

	    set key  [list bond $t1 $t2]
	    set key2 [list bond $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    if {! ([info exists params($key)] ||
		   $t1 == "*" || $t2 == "*")} {
		cgCon -warn "getbondparam undefined parameter $args"
		return -code ok
	    }

	    return -code ok [array get params "$key"]
	}

	angle {

	    if {$narg < 4} {
		cgCon -err "getangleparm missing argument: getangleparam <type1> <type2> <type3>"
		return -code ok
	    }

	    lassign $args command t1 t2 t3

	    set key  [list angle $t1 $t2 $t3]
	    set key2 [list angle $t3 $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    if {! ([info exists params($key)] ||
		   $t1 == "*" || $t2 == "*" || $t3 == "*")} {
		cgCon -warn "getangleparam undefined parameter $args"
		return -code ok
	    }

	    return -code ok [array get params "$key"]
	}

	dihedral {

	    if {$narg < 5} {
		cgCon -err "getdihedralparam missing argument: dihedral t1 t2 t3 t4"
		return -code ok
	    }

	    lassign $args command t1 t2 t3

	    set key  [list dihedral $t1 $t2 $t3 $t4]
	    set key2 [list dihedral $t4 $t3 $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    if {! ([info exists params($key)] ||
		   $t1 == "*" || $t2 == "*" ||
		   $t3 == "*" || $t4 == "*") } {
		cgCon -warn "get_param dihedral undefined parameter $args"
		return -code ok
	    }
	    return -code ok [array get params "$key"]
	}

	pair {

	    if {$narg < 3} {
		cgCon -err "getpairparam missing argument: getpairparam <type1> <type2>"
		return -code ok
	    }

	    lassign $args command t1 t2

	    set key  [list pair $t1 $t2]
	    set key2 [list pair $t2 $t1]

	    ## Check for key existance in reverse order
	    if {[info exists params($key2)]} {set key $key2}

	    if {! ([info exists params($key)] ||
		   $t1 == "*" || $t2 == "*")} {
		cgCon -warn "get_param pair, undefined parameter $args"
		return -code ok
	    }

	    return -code ok [array get params "$key"]
	}

	default {
	    cgCon -err "get*param unknown command: $command. bond, angle, dihedral, pair"
	    return -code error
	}
    }

    return -code ok
}

## Generate pair parameters via selected combination rules
proc CGit::Legacy::param_pair_combine {rule atype1 atype2} {

    lassign [get_param pair $atype1 $atype1] key1 nbparams1
    lassign [get_param pair $atype2 $atype2] key2 nbparams2

    if {$key1 == "" || $key2 == ""} {
	return {UNK -1.0 -1.0}
    }

    lassign $nbparams1 pottype1 epsilon1 sigma1
    lassign $nbparams2 pottype2 epsilon2 sigma2

    ## Check if we have a water interaction
    if {[string match $pottype1 "lj12_4"] ||
	[string match $pottype2 "lj12_4"]} {
	set pottype "lj12_4"
    } else {
	set pottype "lj9_6"
    }

    ## Return the correct interaction type and calculated values
    if {[string match "geometric" $rule]} {
	set s [expr {sqrt(double($sigma1) * double($sigma2))}]
	set e [expr {sqrt(double($epsilon1) * double($epsilon2))}]
	return [list $pottype $e $s]
    } elseif {[string match "arithmetic" $rule]} {
	set s [expr {0.5 * (double($sigma1) + double($sigma2))}]
	set e [expr {0.5 * (double($epsilon1) + double($epsilon2))}]
	return [list $pottype $e $s]
    } elseif {[string match "lb" $rule]} {
	set s [expr {0.5 * (double($sigma1) + double($sigma2))}]
	set e [expr {sqrt(double($epsilon1) * double($epsilon2))}]
	return [list $pottype $e $s]
    } else {
	cgCon -warn "Unknown combination rule"
	return {UNK -1.0 -1.0}
    }
}

## Given a topology array, convert to map format
## Do this transparently when the user reads a legacy top file
proc CGit::Legacy::topo2map {fin {fout {stdout}} {flag 0}} {

    variable sys
    array unset topo *

    ## Read the top file
    read_top $fin topo

    ## Preamble
    set cmd {}
    lappend cmd "if \{!\[namespace exists ::CGit\]\} \{\n"
    lappend cmd "namespace eval ::CGit:: \{\n"
    lappend cmd "variable maptest 1;\n"
    lappend cmd "source common.tcl;\n"
    lappend cmd "\}\n\}\n"

    ## Write out records for each unique residue type

    ## Get data from topo array
    set resname [lsort -unique $topo(resname)]

    foreach r $resname {

	## A proc for each residue type to keep things nice
	lappend cmd "proc ::CGit::map_$r \{\} \{\n"
	lappend cmd "variable map;\n"

	if {$flag == 0} {
	    lappend cmd "puts \"Loading CG mappings from \[info script\]\";\n"
	}

	## Get indices associated with resname of interest
	set index [lsearch -all $topo(resname) $r]

	## Individual bead Properties
	lappend cmd "set map(\[list segname $r\]) \{\};\n"
	set keys {type name charge mass}
	foreach p $keys {
	    lappend cmd "set map(\[list $p $r\]) \{\n"
	    foreach i $index {
		lappend cmd "[lindex $topo($p) $i]\n"
	    }
	    lappend cmd "\};\n"
	}

	## Bond List
	set bnl {}
	foreach b $topo([list bond $r]) {
	    lassign $b idx1 idx2
	    set name1 [lindex $topo(name) $topo($r,$idx1)]
	    set name2 [lindex $topo(name) $topo($r,$idx2)]
	    lappend bnl $name1 $name2
	}

	lappend cmd "set map(\[list bonds $r\]) \{\n"

	foreach {atype1 atype2} $bnl {
	    lappend cmd "\{$atype1 $atype2\}\n"
	}

	lappend cmd "\};\n"

	## Angle List
	set anl {}
	foreach a $topo([list angle $r]) {
	    lassign $a idx1 idx2 idx3
	    set name1 [lindex $topo(name) $topo($r,$idx1)]
	    set name2 [lindex $topo(name) $topo($r,$idx2)]
	    set name3 [lindex $topo(name) $topo($r,$idx3)]
	    lappend anl $name1 $name2 $name3
	}

	lappend cmd "set map(\[list angles $r\]) \{\n"

	if {$anl == {}} {
	    lappend cmd "\{auto\}\n"
	} else {

	    foreach {atype1 atype2 atype3} $anl {
		lappend cmd "\{$atype1 $atype2 $atype3\}\n"
	    }
	}

	lappend cmd "\};\n"
	lappend cmd "return;\n\}\n"
    }

    # Commands to load the proc when the file is sourced
    foreach r $resname {
	lappend cmd "::CGit::map_$r;\n"
	lappend cmd "::CGit::cleanup_map $r;\n"
    }

    lappend cmd "if \{ \$::CGit::maptest \} \{\n"
    lappend cmd "foreach r \{$resname\} \{\n"
    lappend cmd "::CGit::map_stats \$r;\n"
    lappend cmd "::CGit::checkbonds \$r;\n"
    lappend cmd "\}\n\}\n"

    ## Write to stdout/file or load it directly?
    if {$flag} {

	if {[catch {eval [join $cmd]} err]} {
	    cgCon -err "Unable to load map: $err"
	    return -code error
	}

    } else {

	set fid 0
	if {![string equal $fout "stdout"]} {

	    if {[catch {open $fout w} fid]} {
		cgCon -err "mapwrite: Unable to open $fout for writing"
		return -code error $sys(ERROR)
	    }

	    puts "Writing map to $fout"

	} else {
	    set fid stdout
	}

	puts $fid [join $cmd "\n"]

	if {![string equal $fid stdout]} {
	    close $fid
	}
    }

    return -code ok
}

## Write out an importable map for the current selection. Useful if we
## have defined a new molecule or have a PSF file we want to
## convert. We still need to define the atomic makeup of the beads
## afterwards.
proc CGit::Legacy::sel2map {sel {fname {stdout}} {flag 0}} {

    variable sys

    set molid [$sel molid]
    set seltext [$sel text]

    ## Preamble
    set cmd {}
    lappend cmd "if \{!\[namespace exists ::CGit\]\} \{\n"
    lappend cmd "namespace eval ::CGit:: \{\n"
    lappend cmd "variable maptest 1;\n"
    lappend cmd "source common.tcl;\n"
    lappend cmd "\}\n\}\n"

    ## Write out records for each unique residue type in selection

    set resname [lsort -unique -index 0 [$sel get {resname residue}]]
    set resids  [lsearch -all -index 1 -subindices -inline $resname *]
    set resname [lsearch -all -index 0 -subindices -inline $resname *]

    foreach r $resname rid $resids {

        set sel2 [atomselect $molid "($seltext) and residue $rid"]

        if {[$sel2 num] == 0} {continue}

        ## A proc for each residue type to keep things nice
        lappend cmd "proc ::CGit::map_$r \{\} \{\n"
        lappend cmd "variable map;\n"

        if {$flag == 0} {
            lappend cmd "puts \"Loading CG mappings from \[info script\]\";\n"
        }

        ## Individual bead Properties
        set keys {segname type name charge mass}
        foreach p $keys {
            lappend cmd "set map(\[list $p $r\]) \{\n"
            foreach v [$sel2 get $p] {
                lappend cmd "$v\n"
            }
            lappend cmd "\};\n"
        }

        ## Bond List
        set bl [topo -sel $sel2 getbondlist]
        set bnl {}
        set names [$sel get {index name}]
        foreach b $bl {
            lassign $b idx1 idx2
            set name1 [lindex [lsearch -inline -exact -index 0 $names $idx1] 1]
            set name2 [lindex [lsearch -inline -exact -index 0 $names $idx2] 1]
            lappend bnl $name1 $name2
        }

        lappend cmd "set map(\[list bonds $r\]) \{\n"

        foreach {atype1 atype2} $bnl {
            lappend cmd "\{$atype1 $atype2\}\n"
        }

        lappend cmd "\};\n"
        $sel2 delete
    }

    lappend cmd "return;\n\}\n"

    # Commands to load the proc when the file is sourced
    foreach r $resname {
        lappend cmd "::CGit::map_$r;\n"
        lappend cmd "::CGit::cleanup_map $r;\n"
    }

    lappend cmd "if \{ \$::CGit::maptest \} \{\n"
    lappend cmd "foreach r \{$resname\} \{\n"
    lappend cmd "::CGit::map_stats \$r;\n"
    lappend cmd "::CGit::checkbonds \$r;\n"
    lappend cmd "\}\n\}\n"

    ## Write to stdout/file or load it directly?
    if {$flag} {

        if {[catch {eval [join $cmd]} err]} {
            cgCon -err "Unable to load map: $err"
            return -code error
        }

    } else {

        set fid 0
        if {![string equal $fname "stdout"]} {

            if {[catch {open $fname w} fid]} {
                cgCon -err "mapwrite: Unable to open $fname for writing"
                return -code error $sys(ERROR)
            }

            puts "Writing map to $fname"

        } else {
            set fid stdout
        }

        puts $fid [join $cmd "\n"]

        if {![string equal $fid stdout]} {
            close $fid
        }
    }

    return -code ok
}

## Convert PSF file to map file
## FIX ME
proc ::CGit::Legacy::psf2map {psffile {mapfile stdout}} {

    variable sys
    array set map {}

    ## Open and parse the psf file
    if {[catch {open $psffile r} fid]} {
        cgCon -err "mapwrite: Unable to open $fname"
        return -code error $sys(ERROR)
    }

    set data [read $fid]
    close $fid

    ## Search for "!NBOND" identifier
    set idx [lsearch -glob $data "*!NBOND*"]
    if {$idx == "-1"} {
        puts stderr "Missing \"!NBOND\" identifier, malformed psf file?"
        return -code error
    }

    #Get number of bonds
    lassign [lindex $data $idx-1] nbonds

    ## Read the next nbonds
    set bonds {}
    set bl [lrange $data $idx [expr {$idx + $nbonds * 2 + 1}]]
    foreach {id1 id2} $bl {
        lappend bonds [list $id1 $id2]
    }

    set data [split $data "\n"]

    ## search for "!NATOM" identifier
    set idx [lsearch -glob $data "*!NATOM*"]
    if {$idx == "-1"} {
        puts stderr "Missing \"!NATOM\" identifier, malformed psf file?"
        return -code error
    }

    #Get number of atoms
    lassign [lindex $data $idx] natoms

    ## Read the next natom lines
    for {set i 1} {$i <= $natoms} {incr i} {
        lassign [lindex $data $idx+$i] atomid segname resid resname\
            atomname atomtype charge mass

        lappend m_id($resname)       $atomid
        lappend m_name($resname)     $atomname
        lappend m_type($resname)     $atomtype
        lappend m_charge($resname)   $charge
        lappend m_mass($resname)     $mass
    }

    ## Preamble
    set cmd {}
    lappend cmd "if \{!\[namespace exists ::CGit\]\} \{\n"
    lappend cmd "namespace eval ::CGit:: \{\n"
    lappend cmd "variable maptest 1;\n"
    lappend cmd "source common.tcl;\n"
    lappend cmd "\}\n\}\n"

    ## Get unique resnames
    set resnames [lsort -unique [array names m_atomid *]]

    ## Foreach unique resname, atomname generate the map
    set idx {}
    foreach r $resname {
        set names [lsort -unique $m_name($r)]
        foreach n $names {
            lappend idx [lsearch $m_name($r) $n]
        }

        ## A proc for each residue type to keep things nice
        lappend cmd "proc ::CGit::map_$r \{\} \{\n"
        lappend cmd "variable map;\n"
        lappend cmd "puts \"Loading CG mappings from \[info script\]\";\n"

        foreach p {name type charge mass} {
            lappend cmd "set map(\[list $p $r\]) \{\n"
            foreach v $idx {
                lappend cmd "m_$p($v)\n"
            }
            lappend cmd "\};\n"
        }

        lappend cmd "set map(\[list bonds $r\]) \{\n"

        return -code ok

    }
}

# +------------------------------------------------------------+
# | Guess the relative weights for the COM calculations        |
# | If an atom is listed twice in a bead, it is assumed        |
# | to be "split" between beads, and thus the relative weights |
# | for the atom in each bead would be 0.5.                    |
# +------------------------------------------------------------+

proc ::CGit::Legacy::mapWeights {resname} {
    variable map

    foreach r $resname {
        ## Check if the user has requested to guess the weights
        if {[catch {dict get $map $r weights}] ||
            [lsearch [dict get $map $r weights] "auto"] >= 0} {

            set weights {}

            ## Get a list of all the unique atoms
            if {[catch {dict get $map $r map} atoms]} {
                cgCon -error "Missing map for residue $r"
                return -code error
            }

            set atoms [join $atoms]
            set uniq_atoms [lsort -unique $atoms]

            set n_atoms [llength $atoms]
            set n_uniq [llength $uniq_atoms]

            ## if the number of atoms equals the number of unique
            ## atoms, we have no duplicates, all weights are 1.0
            if {$n_atoms == $n_uniq} {
                set weights [lrepeat $n_atoms 1.0]
            } else {
                ## Tabulate the weights for each unique atom
                foreach a $uniq_atoms {
                    set n [llength [lsearch -all $atoms $a]]
                    lappend weights [expr 1.0/$n]
                }
            }

            ## Get a list of all the beads
            if {[catch {dict get $map $r map} beads]} {
                cgCon -error "Missing map for residue $r"
                return -code error
            }

            ## Construct the weight list
            set ww {}
            foreach b $beads {
                set w {}
                foreach a $b {
                    set idx [lsearch $uniq_atoms $a]
                    lappend w [lindex $weights $idx]
                }
                lappend ww $w
            }

            ## Add it to the dictionary
            dict set map $r weights $ww

        } else {

            ## make sure we have a weight for each atom
            ## if we're not guessing.

            ## Get a list of all the atoms
            if {[catch {dict get $map $r map} atoms]} {
                cgCon -error "Missing map for residue $r"
                return -code error
            }

            ## Get a list of all the weights
            if {[catch {dict get $map $r weights} weights]} {
                cgCon -error "Missing weights for residue $r"
                return -code error
            }

            foreach a $atoms w $weights {
                if {[llength $a] != [llength $w]} {
                    cgCon -error "Missing or too many weighting coefficents for residue $r"
                    return -code error
                }
            }
        }
    }

    return -code ok
}
