# Testing initial value computations with auxiliary variables (2 period lags, diff),
# nonlinear functions and an initval file 
var x, y;
varexo e;

parameters a1, rho1, rho2;
rho1 =  0.8;
rho2 = -0.2;
a1 = 0.11;

model;
  y = rho1*y(-1) + rho2*y(-2) + e;
  x = diff(y*y + a1*y^2 - exp(y) + exp(ln(y)) + sin(y) - cos(y) + y/y(-1));
end;

initval_file(datafile='data_nl.csv');
 
