#!/usr/bin/tclsh8.5

# Make a coarse-grained trajectory from an all-atom trajectory.

proc trajectory {} {

  ## Load the cg-tools package  
  package require cg

  ## Load Parameters
  cg readparam ../common/par_all1_sdk.json 
  cg readtop   ../common/top_all1_sdk.json 

  ## Load the all-atom trajectoty
  set molid [mol new ./water-heptane.psf]
  mol addfile ./water-heptane.xtc waitfor all

  ## Coarse-grain the trajectory
  set cgid [cg -molid $molid maptraj -all]

  ## Set the bond/angle topologies
  ## apply the bead properties
  cg -molid $cgid setbonds
  cg -molid $cgid setangles
  cg -molid $cgid reanalyzemol

  ## Write out the coarse-grained trajectory
  animate write psf water-heptane.cg.psf $cgid
  animate write dcd water-heptane.cg.dcd $cgid

  ## write out the lammps inputs
  cg -molid $cgid writelammpsdata water-heptane.data full
  cg -molid $cgid writelammpsparam water-heptane.param both
}

if {1} {trajectory}
