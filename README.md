# Simulation Codes in PhD Project.
## Thesis Topic: Fracture dynamics of correlated percolation on polymer networks.

My PhD project is to study fracture (crack) dynamics on a polymer network. The practical goal is to predict the lifetime of the polymer network, but as it can be modelled on a random network, this is basically a mathematical random graph problem that can be used to predict the lifetime of a given network. As the random fracture process is a Markov process, the dynamical process is numerically simulated in Python and implemented with a kinetic Monte Carlo method.

During the dynamics, we focused on locating a spanning or percolating clusters in the network. This is solved by implementing an efficient tree-based algorithm. 

This code is a whole process of how the stochastic dynamics of the network is simulated in Python.

##
SL_gamma_redistribution_bond_UpDown_or_LeftRight.py is the main script that simulates the dynamic process.

