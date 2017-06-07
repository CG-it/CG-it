# CG-it

CG-it is a VMD package to simplify creating coarse grained (CG) topologies for molecular dynamics simulations.  It currently supports CG force fields that follow the "SDK" parameterization procedure from Shinoda et al, Mol Sim, 33:27-36 (2007)

Developed by Chris MacDermaid <chris.macdermaid@gmail.com> and Giacomo Fiorin <giacomo.fiorin@gmail.com>

Special thanks to Axel Kohlmeyer for advice and writing [TopoTools](https://github.com/akohlmey/topotools)

## Install
1. From the main package directory, issue `make install` to compile and copy the relevant files to `${HOME}/.vmdplugins/` folder.
2. Add the following to your `${HOME}/.vmdrc` file:
```
set auto_path [linsert $auto_path 0 [file join $env(HOME) .vmdplugins]]
```
3. Start VMD and issue the command (or add it to `${HOME}/.vmdrc` for convenience):
```
package require cg
```


## Example input

To read the parameters database file (available [here](https://github.com/CG-it/ffdb-sdk)):
```
cg readparam parameters.json
```
To read the topology database file (available [here](https://github.com/CG-it/ffdb-sdk)):
```
cg readtop topology.json
```
To convert the current VMD molecule from atomistic to CG:
```
cg map
```
To assign the types of bonded interactions (not all commands may be needed):
```
cg setbonds
cg setangles
cg setdihedrals
cg setimpropers
```
To add water to the system if not already present (unit cell must be defined):
```
cg solvate
```
To finalize the CG topology:
```
cg reanalyzemol
```
To get the full help invoke the `cg` command without arguments:
```
cg help
```

