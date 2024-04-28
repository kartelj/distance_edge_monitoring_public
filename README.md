# Instructions for use

The repository contains all the code, binaries, instances and run scripts required to reproduce the results reported in the paper.
There are three folders:

1. code -- containing the IP implementation that depends on the CPLEX solver library. We used the IBM ILOG CPLEX Studio 221 version of the library, and it should probably work with the later versions as well.
2. instances -- with all tested instances divided into two subfolders: instances/random and instances/regular.
3. test -- with the script run_cplex.cmd, which is used to call the cplex binary (there are already prebuilt Windows binaries in test folder: dem_cplex.exe and dem_cplex.pdb) and reproduces the reported results.
