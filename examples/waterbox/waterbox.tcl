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

## CG and replicate a waterbox and make lammps inputs

proc waterbox {} {

    package require cg

    ## Load the molecule
    set molid [mol new ./waterbox.eq300.pdb]

    ## Load parameters for water
    cg readparam ../common/par_all1_sdk.json
    cg readtop   ../common/top_all1_sdk.json

    ## By default, cg works on the top molecule and
    ## all atoms within that molecule. Use the -sel
    ## or -molid flags to specify particular atoms/mols

    ## make CG model
    set newmol [cg map]

    ## Assign bead properties
    cg reanalyzemol

    ## Set the approximate box dimensions
    set sel [atomselect $newmol "all"]
    set box [vecinvert [vecsub {*}[measure minmax $sel]]]
    molinfo $newmol set {a b c} $box
    $sel delete

    ## Replicate the box 2x2x2
    set newmol [::TopoTools::replicatemol $newmol 2 2 2]

    ## Write out lammps files
    cg -molid $newmol writelammpsdata waterbox-cg.data full
    cg -molid $newmol writelammpsparam waterbox-cg.param cpu
}

if {1} {waterbox}
