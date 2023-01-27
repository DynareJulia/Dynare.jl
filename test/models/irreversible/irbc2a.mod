// removing the constraint

var lk, lc, i, A, a, mu;
varexo epsilon;
parameters beta, alpha, delta, rho, Astar, sigma, L;

beta    =  0.990;
alpha   =  0.450;
delta   =  0.020;
rho     =  0.980;
Astar   =  1.000;
sigma   =  1.000;
L       =  1.000;

model;
  a = rho*a(-1) + epsilon;
  A = Astar*exp(a);
  exp(lc*(-sigma)) - mu = beta*(exp(lc(+1)*(-sigma))*(alpha*A(+1)*exp(lk*(alpha - 1))*L^(1 - alpha)
                    + 1 - delta) - mu(+1)*(1 - delta));
  exp(lk) + exp(lc) = A*exp(lk(-1)*alpha)*L^(1 - alpha) + (1 - delta)*exp(lk(-1));
  i = exp(lk)-(1 - delta)*exp(lk(-1));
[mcp = 'i > 0' ]
  mu = 0;

end;

steady_state_model;
  a = 0;
  mu = 0;
  A = Astar;

  lk = log(((1/beta - 1 + delta)/(alpha*A*L^(1 - alpha)))^(1/(alpha - 1)));
  i = delta*exp(lk);
  lc = log(A*exp(lk*alpha)*L^(1 - alpha) - i);
end;

shocks;
 var epsilon;
 periods 10;
 values -2;
end;

steady;

perfect_foresight_setup(periods=400);
perfect_foresight_solver(lmmcp);

