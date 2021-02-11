# Simulation Codes in PhD Project.
## Thesis Topic: Fracture dynamics of correlated percolation on polymer networks.

My PhD project targeted studying fracture (crack) dynamics on a polymer network. A reader can find my article (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.101.042603) and my thesis in this repository to have a general understanding.

The practical goal is to predict the lifetime of the polymer network. But as the fractures on a polymer network is a random network problem, the physical dynamics can be mapped to a mathematical random graph.


As the random fracture process is a Markov process also, in which the random fracture events exhibit correlations. Thus, the dynamical process is numerically simulated with a kinetic Monte Carlo method in Python.

During the dynamics, we focused on locating a spanning (or percolating) cluster in the network. This is again a math problem called percolation theory. In a numerical simulation, we can let the computer or the program identify such percolating cluster by an efficient tree-based algorithm by Newman (https://link.aps.org/doi/10.1103/PhysRevE.64.016706). 

This program was created by me in Python to simulate a whole random dynamic process. But remember, as a one probable process is only one possible copy. I have to simulate many copies for each set of parameters.

33

#### What is fracture dynamics on polymer networks? How do I map the dynamics of the polymer networks onto a random graph problem? What is a percolating cluster? Look at the below image!
![image](https://github.com/yulew/PhDProject/blob/main/imgs/Mapping.png)

##
SL_gamma_redistribution_bond_UpDown_or_LeftRight.py is the main script that simulates the dynamic process.

