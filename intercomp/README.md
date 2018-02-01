# README #

InterComp is a software for the sequence-order alignment and superposition of biological molecules.

Currently InterComp supports two input formats:

* PDB
* mol2

## Installation ##

To compile InterComp just clone this repo and run the command:

```
make
```

A binary named InterComp will be generated in the cloned folder.

## Usage examples ##

When aligning two PDB files with default optimization parameters, InterComp is called with:

```
./InterComp -PDB path/to/file1.pdb path/to/file2.pdb
```

The order of the two files following the "-PDB" flag is irrelevant.

In this case, InterComp only outputs the structural and sequence score of the alignment.

If a superposition of the two molecules is needed, use:

```
./InterComp -PDB path/to/file1.pdb path/to/file2.pdb -super
```

This will generate two PDB files:

* "molA.pdb", which contains the CA atoms from the smallest of the two input molecules
* "molB.pdb", which contains only those CA atoms from the largest of the input molecules that superimpose to the atoms in molA.pdb. These are also rotated and translated accordingly onto the atoms in "molA.pdb".
