# Pair potential pair_lj_exp_att for LAMMPS (2016) version (included)

This repository contains a custom pair potential, that implements an adjustable Weeks-Chandler-Anderson potential (lj_cut) in addition to a short ranged attractive tail based upon an exponential decay, the Lammps 2016 stable release to which it is fully compatiable with, and a simple working example (2 beads in a box).

### Installing

First untar the lammps-stable-2016.tar.gz at the terminal :
```
tar -xvzf lammps-stable-2016.tar.gz
```
which unpacks the LAMMPS files into a directory ***lammps-5Nov16/***. Then execute these commands:
```
cp pair_lj* lammps-5Nov16/src/
```
then
```
cd lammps-5Nov16/src/
make clean-all
make purge
make serial
```

### Test example

This might take a few moments. Then once that has done successfully i.e. no errors or warnings appear. Execute:
```
cp ../../config.txt $PWD
cp ../../in.input_script $PWD
```
Then to check everything works we run a test simulation with 2 beads. You simply run this on the command line:
```
./lmp_serial -in in.input_script
```
If the last line is `=== pair...  works ===` then you are ready to use pair_lj_exp_att.

### Usage

The `pair_style`:

```
pair_style lj/exp/att <cut> <dia> <epp> <WCA-switch>
```
where <cut> is the cutoff for the total interaction i.e. any beads at or beyond this distance with have zero potential energy, <dia> is the hard sphere diameter (specifically the cutoff for the Weeks-Chandler-Anderson potential) and is smaller than <cut>, <epp> is the minimum of the pair potential (at <dia), and <WCA-switch> switches the Weeks-Chandler-Anderson potential on/off, e.g., if you wanted to use another repulsive energy potential instead.

The `pair_coeff`:
```
pair_coeff <type_1> <type_2> <eps> <lambda>
```
where <type_1> and <type_2> are particle types in LAMMPS, <eps> is the strength of the Weeks-Chandler-Anderson potential, and <lambda> is the decay length (normally set to the particle diameter).

# ljexpatt-pairpotential
