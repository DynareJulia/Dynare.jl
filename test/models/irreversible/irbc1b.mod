// removing the constraint

var ek, ec, i, A, a, mu;
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
  exp(ec)^(-sigma) - mu = beta*(exp(ec(+1))^(-sigma)*(alpha*A(+1)*exp(ek)^(alpha - 1)*L^(1 - alpha)
                    + 1 - delta) - mu(+1)*(1 - delta));
  exp(ek) + exp(ec) = A*exp(ek(-1))^alpha*L^(1 - alpha) + (1 - delta)*exp(ek(-1));
  i = exp(ek)-(1 - delta)*exp(ek(-1));
[mcp = 'i > 0' ]
  mu = 0;

end;

steady_state_model;
  a = 0;
  mu = 0;
  A = Astar;

  ek = log(((1/beta - 1 + delta)/(alpha*A*L^(1 - alpha)))^(1/(alpha - 1)));
  i = delta*exp(ek);
  ec = log(A*exp(ek)^alpha*L^(1 - alpha) - i);
end;

shocks;
 var epsilon;
 periods 10;
 values -1;
end;

steady;

perfect_foresight_setup(periods=400);
perfect_foresight_solver(lmmcp);

