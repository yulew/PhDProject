# Simulation Codes in PhD Project.
## Thesis Topic: Fracture dynamics of correlated percolation on polymer networks.

My PhD project targeted studying fracture (crack) dynamics on a polymer network. A reader can read my article (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.101.042603) and my thesis in the repository "article and PhDthesis".

The practical goal of the project is to forcast the lifetime of polymer networks. But as the fractures on a polymer network is a random network problem, the physical dynamics can be mapped to a mathematical random graph.


As the random fracture process is a Markov process also, in which the random fracture events exhibit correlations in the time-series processes. Thus, the fracture dynamics were numerically simulated with a kinetic Monte Carlo method in Python.

During the dynamics, we focused on locating a spanning (or percolating) cluster in the network. This is again a math problem called percolation theory. In a numerical simulation, we can let the computer or the program identify such percolating cluster by an efficient tree-based algorithm by Newman (https://link.aps.org/doi/10.1103/PhysRevE.64.016706). 

I created this program in Python to simulate a whole random dynamic process. But remember, as a one probable process is only one possible copy. I have to simulate many copies for each set of parameters.

##

What is fracture dynamics on polymer networks? How do I map the dynamics of the polymer networks onto a random graph problem? What is a percolating cluster? Check out the below image!

![image](https://github.com/yulew/PhDProject/blob/main/imgs/Maping.png)

##
What is a percolation theory? It is a theory that you want to find the probability of generating a spanning/percolating cluster. Like this!
<div align=center><img width="550" height="200" src="https://github.com/yulew/PhDProject/blob/main/imgs/percolation.png">


##
What is a (rejection-free) kinetic Monte Carlo simulation method? It is a random sampling method, but the next event is dependent on its current state. The below describe the simulation process of the fracture dynamics!
<div align=center><img width="550" height="500" src="https://github.com/yulew/PhDProject/blob/main/imgs/Monte_Carlo.png">

##
Check out the animations of the dynamic process of the polymer network!

An example of a low correlated fracture dynamics system, more randomly: https://github.com/yulew/PhDProject/blob/main/Animations/animation_perco_low_corr.mp4

An example of a highly correlated dynamics system, the large "accumulating" crack starts so early: 
https://github.com/yulew/PhDProject/blob/main/Animations/animation_high_corr.mp4

But you should know, even for a completely same initial condition, the random dynamic processes can be completely different.


##

SL_gamma_redistribution_bond_UpDown_or_LeftRight.py is the main script that simulates the dynamic process.

