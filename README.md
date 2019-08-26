# Optimal transportation networks via a phase field approximation of the generalized urban planning problem 

## General information 

The code implements the method proposed in [[1]](https://arxiv.org/abs/1805.11399), for solving the generalized urban planning problem [[2]](http://www.numdam.org/article/COCV_2005__11_1_88_0.pdf) via a phase field approximation approach. 

The problem was implemented using [Matlab](https://www.mathworks.com/products/matlab.html), Version R2018b. 


## Run example 

To run the code, open Matlab and execute Example.m. This will compute a simple example of transport from one point source to two point sinks with equal mass. For an other example or a different setting, the following parameters can be changed within the example file:

	n: Image size (n by n grid) 
	savename: Savename 
	epsilonStart: Initial epsilon 
	epsilonEnd: Final epsilon 
	maxIter: Maximal number of iterations
	tol: Iteration error tolerance  
	a0: Parameter a0
	a: Parameter a 
	b: Parameter b 
	example: Example (see getExample.m)
	initSol: Initial solution given (true or false)
	savenameInit: Filename of initial solution (filename/savename of given initial solution)


[1] Luca Ferrari, Carolin Rossmanith, and Benedikt Wirth. Phase field approximations of branched transportation problems. Preprint, 2018. 

[2] Alessio Brancolini and Giuseppe Buttazzo. Optimal networks for mass transportation problems. *ESAIM Control Optim. Calc. Var.*, 11(1):88-101 (electronic), 2005.



