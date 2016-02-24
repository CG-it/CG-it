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

## Number of unique resdiues in selection
proc CGit::nres {sel} {return [llength [lsort -unique [$sel get residue]]]}
