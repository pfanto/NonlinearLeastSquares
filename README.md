# NonlinearLeastSquares
Nonlinear least-squares fitting solvers

This code minimizes the loss function L = || y - f(x|a) ||^2 for the parameters a of a nonlinear function f.
The main solver uses the Levenberg-Marquardt algorithm (see Press et al., Numerical Recipes).
The other solver uses the standard linearization method with regularization.





