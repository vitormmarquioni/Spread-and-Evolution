READ ME - EpidemicEvolution.f90

////////////////////////////////////////////////////////

This program evaluates the dynamics of an epidemic on a network
and calculates the pairwise distance between all the viruses,
those that are still active and those that are no longer active.

The code is written in fortran 90.
I use to write a Wolfram Language notebook to write the input
files and to analize the outputs

The details of the epidemic spread can be found in the paper
by Marquioni and de Aguiar, 2020, "Quantifying the effects
of quarantine using an IBM SEIR model on scalefree networks",
https://doi.org/10.1016/j.chaos.2020.109999

In this version, the transmission probability is kept constant,
but it can be changed as function of time, as it has been done
in the mentioned paper.

Any comments or suggestions can be sent to by email:
vimarqmon@gmail.com

////////////////////////////////////////////////////////
This program needs the following input files:
(In the repository you can find example input files)
-------------------Inputs------------------

---> "auxi.in"
-------------------------------------------
- include the parameters in the following order:
1)population size (integer)
2)genome size times 2 (integer)
3)time horizon (integer)
4)transmission probability (real)
5)recovery probability (real)
6)mutation rate (real)
7)initial infected individual (integer)
8)name the network file (string)

---> "auxi2.in"
-------------------------------------------
- this files says the number of modules the network has and
each size. Of course, a network with only one component can
be described in this code as a single module network.
- this code allows a maximum of 10 modules, but this can be
easily changed if one changes the names of the outputs.
- The file is written in the following order:
1)number of modules
2)size of module 1
3)size of module 2
...
////////////////////////////////////////////////////////
---> "seed.in"
-------------------------------------------
- this is a file with 12 integer numbers separated by blank
spaces.
-this file is updated at the end of each run.
////////////////////////////////////////////////////////
---> "incub.in"
-------------------------------------------
- this file has a list of the incubation periods of every
individual, i.e., if an individual becomes exposed, it will be
exposed for this time before becoming infected. Once the
inbution time distribution is not simple, it is easier to
generate these times out of the main program and then import
them.
- the main program reads a vertical lists, where each time is
in a different line.
////////////////////////////////////////////////////////
---> adjacency "network" file
-------------------------------------------
- the program reads de adjacency matrix of the contact network.
- this file, whose name is given in "auxi.in", is written has.
- the adjacency matrix must be given in a way that the individuals
are numbered according to its module.
- this program considers undirected graphs.
- The file is written in the following format:
a(1,2)
a(1,3)
...
a(1,n)
a(2,3)
a(2,4)
...
a(2,n)
...
a(n-1,n)
////////////////////////////////////////////////////////

This program generates the following output files:
------------------Outputs------------------

---> "OUT4.txt"
-------------------------------------------
- information file about the final state of the simulation.
////////////////////////////////////////////////////////
---> "ancestor.txt":
-------------------------------------------
- information about who infected who and when it happened
- written in the format:
(individual that was infected),(who infected),(when was infected)
////////////////////////////////////////////////////////
---> "extinction.txt"
-------------------------------------------
- information about when an individual recovered
- written in the format:
(when),(who)
////////////////////////////////////////////////////////
---> "phylo.out"
-------------------------------------------
- infection tree
- written in the format:
(infecting individual)->(infected individual)
////////////////////////////////////////////////////////
---> "fort.1000+i"
-------------------------------------------
- distance matrix among viruses at time i
- these files are written in the following format:
dist(1,2)
dist(1,3)
...
dist(1,n)
dist(2,3)
dist(2,4)
...
dist(2,n)
...
dist(n-1,n)
////////////////////////////////////////////////////////
---> "fort.100+i"
-------------------------------------------
- state vector at time i
////////////////////////////////////////////////////////
Each line j in every file below has a number corresponding
to time step j-1.
-------------------------------------------
---> "fort.11": total number of susceptibles (%).
---> "fort.12": total number of infecteds+exposeds (%).
---> "fort.13": total number of recovereds (%).
-------------------------------------------
---> "susModK.out": number of susceptibles in module K (%).
---> "infModK.out": number of infecteds+exposeds in module K (%).
---> "recModK.out": number of recovered in module K (%).
-------------------------------------------
