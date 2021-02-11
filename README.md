# Simulation Codes in PhD Project.
## Thesis Topic: Fracture dynamics of correlated percolation on polymer networks.

My PhD project targeted studying fracture (crack) dynamics on a polymer network. A reader can find my article (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.101.042603) and my thesis in this repository to have a general understanding.

The practical goal is to predict the lifetime of the polymer network. But as the fractures on a polymer network is a random network problem, the physical dynamics can be mapped to a mathematical random graph.


As the random fracture process is a Markov process also, in which the random fracture events exhibit correlations. Thus, the dynamical process is numerically simulated with a kinetic Monte Carlo method in Python.

During the dynamics, we focused on locating a spanning (or percolating) cluster in the network. This is again a math problem called percolation theory. In a numerical simulation, we can let the computer or the program identify such percolating cluster by an efficient tree-based algorithm by Newman (https://link.aps.org/doi/10.1103/PhysRevE.64.016706). 

I created this program in Python to simulate a whole random dynamic process. But remember, as a one probable process is only one possible copy. I have to simulate many copies for each set of parameters.

##

What is fracture dynamics on polymer networks? How do I map the dynamics of the polymer networks onto a random graph problem? What is a percolating cluster? Check out the below image!
![image](https://github.com/yulew/PhDProject/blob/main/imgs/Mapping.png)

##
What is a percolation theory? It is a theory that you want to find the probability of generating a spanning/percolating cluster. Like this!
<div align=center><img width="550" height="200" src="https://github.com/yulew/PhDProject/blob/main/imgs/percolation.png">

##
What is a (rejection-free) kinetic Monte Carlo simulation method? It is a random sampling method, but the next event is dependent on its current state. The below describe the simulation process of the fracture dynamics!
<div align=center><img width="550" height="550" src="https://github.com/yulew/PhDProject/blob/main/imgs/Monte_Carlo.png">

##
Check out the animation of the dynamic process!

A low correlated fracture system:

[![Alt text for your video](doc/screenshot_youtube.PNG)](https://www.youtube.com/watch?v=5EEM__L612Y "Put hover text here!")

##

SL_gamma_redistribution_bond_UpDown_or_LeftRight.py is the main script that simulates the dynamic process.

