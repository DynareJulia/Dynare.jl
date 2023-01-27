// removing the constraint

var k, c, i, A, a, mu;
varexo L, epsilon;
parameters beta, alpha, delta, rho, Astar, sigma;

beta    =  0.990;
alpha   =  0.450;
delta   =  0.020;
rho     =  0.980;
Astar   =  1.000;
sigma   =  1.000;

model;
  a = rho*a(-1) + epsilon;
  A = Astar*exp(a);
  c^(-sigma) - mu = beta*(c(+1)^(-sigma)*(alpha*A(+1)*k^(alpha - 1)*L(+1)^(1 - alpha)
                    + 1 - delta) - mu(+1)*(1 - delta));
  k + c = A*k(-1)^alpha*L^(1 - alpha) + (1 - delta)*k(-1);
  i = k-(1 - delta)*k(-1);

//  [ mcp = 'i > 0' ]
  mu = 0;
end;

steady_state_model;
  a = 0;
  mu = 0;
  A = Astar;

  k = ((1/beta - 1 + delta)/(alpha*A*L^(1 - alpha)))^(1/(alpha - 1));
  i = delta*k;
  c = A*k^alpha*L^(1 - alpha) - i;
end;

shocks;
 var epsilon;
 periods 10;
 values -1.0;
end;

initval;
  L = 1;
end;

steady;

perfect_foresight_setup(periods=400);
perfect_foresight_solver;

