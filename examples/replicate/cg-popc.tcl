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

## CG a popc bilayer-water system and make lammps outputs

proc cg_popc {} {

    package require cg
    package require topotools

    ## Load the parameters
    cg readparam ../common/par_all1_sdk.json 
    cg readtop   ../common/top_all1_sdk.json

    ## Load the molecule
    set molid [mol new ./sb.pdb]

    ## By default, cg works on the top molecule and
    ## all atoms within that molecule. Use the -sel
    ## or -molid flags to specify particular atoms/mols

    ## make CG model
    set newmol [cg map]

    ## Assign topology information
    cg setbonds
    cg setangles
    cg reanalyzemol

    ## Set the approximate box dimensions
    set sel [atomselect $newmol "all"]
    set box [vecinvert [vecsub {*}[measure minmax $sel -withradii]]]
    molinfo $newmol set {a b c} $box
    $sel delete

    ## Replicate the system
    set newmol [::TopoTools::replicatemol $newmol 20 20 1]

    ## Write out a pdb and psf for easy work with VMD 
    animate write psf lb-cg.psf
    animate write pdb lb-cg.pdb

    ## Write out lammps files
    cg writelammpsdata lb-cg.data full
    cg writelammpsparam lb-cg.param cpu
}

after idle {
    if {1} {cg_popc}
}
