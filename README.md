# NonlinearLeastSquares
Nonlinear least-squares fitting solvers

This code minimizes the loss function L = || y - f(x|a) ||^2 for the parameters a of a nonlinear function f.
The main solver uses the Levenberg-Marquardt algorithm (see Press et al., Numerical Recipes).
The other solver uses the standard linearization method with regularization.

The code assumes that the function takes as arguments the input x, the number of model parameters npars, and a double pointer array "a" containing all the model parameters.  The number of adjusted parameters in the least-squares fit is m <= npars.  The parameters to be adjusted should be the first m values in the pointer array a.  

A test case of fitting a simple exponential function is included in the source directory.
To perform this test, run the shell commands
$ make
$ ./test

Direct any questions to paul.fanto@yale.edu.


