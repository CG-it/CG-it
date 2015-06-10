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

## Solvate a POPC bilayer

proc solvate {} {

    package require cg

    ## Load the popc bilayer
    set molid [mol new ./popc-cg.psf]
    mol addfile ./popc-cg.pdb molid $molid

    ## Load parameters
    cg readparam ../common/par_all1_sdk.json
    cg readtop ../common/top_all1_sdk.json

    ## Apply bead properties
    cg -molid $molid reanalyzemol

    ## Get the bilayer dimensions
    set sel [atomselect $molid "all"]
    set dim [vecinvert [vecsub {*}[measure minmax $sel]]]

    ## Add an extra 10 angstroms of padding along z-dim on each side of bilayer
    set dim [vecadd $dim {0.0 0.0 20.0}]

    ## Set the box size
    molinfo $molid set {a b c} $dim

    ## Solvate the bilayer (remove waters within 4.0 angstoms of solute)
    set newmol [cg solvate 4.0]

    ## Write it out
    animate write psf popc-solvate.cg.psf
    animate write pdb popc-solvate.cg.pdb

    ## Make lammps input files
    cg -molid $newmol writelammpsdata popc-solvate.data full
    cg -molid $newmol writelammpsparam popc-solvate.param both 
    
    # Cleanup
    $sel delete 
}

if {1} {solvate}
