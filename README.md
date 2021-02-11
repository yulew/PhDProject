# Simulation Codes in PhD Project.
## Thesis Topic: Fracture dynamics of correlated percolation on polymer networks.

My PhD project targeted studying fracture (crack) dynamics on a polymer network. A reader can find my article (DOI: 10.1103/PhysRevE.101.042603) and my thesis in this repository.

The practical goal was to predict the lifetime of the polymer network, but as it can be modelled on a random network, this is basically a mathematical random graph problem that can be used to predict the lifetime of a given network. As the random fracture process is a Markov process, the dynamical process is numerically simulated in Python and implemented with a kinetic Monte Carlo method.

During the dynamics, we focused on locating a spanning (or percolating) cluster in the network. This can be solved by implementing a tree-based algorithm by Newman for efficiency (https://link.aps.org/doi/10.1103/PhysRevE.64.016706). 

This program was created by me in Python to simulate a whole random dynamic process. But remember, as a one probable process is only one possible copy. I have to simulate many copies for each set of parameters.

##
SL_gamma_redistribution_bond_UpDown_or_LeftRight.py is the main script that simulates the dynamic process.

