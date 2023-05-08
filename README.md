# MDPAT

Molecular Dynamics Post-processing and Analysis Toolkit by Michael Jacobs

A high-performance, high-throughput application meant for supercomputers (e.g., Summit and Frontier at the Oak Ridge Leadership Computing Facility).

Designed to be modular with swappable interfaces for computer architecture, trajectory file formats, etc.

## TODO

 - Fix type warnings
 - Update readTrajectories tests
 - Finish serial msd tests
 - Do parallel msd tests
 - Do readInput tests
 - Finish readInput
 - Upgrade permuteDimsParallel for efficiency

## `src/initNode.cpp`

This initializes MPI and the GPUs on the system. Currently, only Summit's architecture is supported, but we will add support for Frontier as well.

## `src/readInput.cpp`

This reads a user file via `stdin` with several commands defining the dump files, atom types, timestep, degree of polymerization, etc.

## `src/readTrajectories.cpp`

This is meant to be a universal tool for reading a LAMMPS dump file with any number of columns representing any values. Additionally, this should allow for the processing of LAMMPS binary dump files for faster processing.

## `src/permuteDims.cpp`

A convenience function meant to permute the dimensions of an array to group elements that should be processed together. For example, in the computation of the mean-squared displacement, the smallest grouping would be all trajectories of a single component of a single atom.

## `src/selectFromArray.cpp`

Not entirely sure if I'll need this in the future... 

## `src/splitValues.cpp`

A convenience function to efficiently split a number of values among the participating processors. This can be used more than once in a script, for example, one could call this over dump files to parallelize reading them, then call it to split up the trajectories by atom to analyze as mentioned in the `permuteDims.cpp` section.

## `src/output.cpp`

Functions to write out data to files, either by columns or as a table.

## `bcastContainers.cpp`

Functions to simplify broadcasting strings and vectors, particularly when reading input.

## `mpiType.hpp`

This was taken from [github](https://gist.github.com/2b-t/50d85115db8b12ed263f8231abf07fa2) to add flexibility to the types used in C-style MPI calls.

## `tools/scatteringVectors.c`

A group of functions to generate a file for describing a fixed collection of scattering vectors over which some scattering analysis should take place.

