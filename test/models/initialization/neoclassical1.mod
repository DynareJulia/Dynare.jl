// convergence to equilibrium

// endogenous variables declaration
var c k;
// exogenous variables declaration
varexo A;

// parameters declaration
parameters alpha beta delta sigma;
alpha = 0.3;
beta = 0.98;
delta = 0.025;
sigma = 1;

// model equations
model;
c^(-sigma) = beta*c(+1)^(-sigma)
               *(alpha*A(+1)*k^(alpha-1)+1-delta);
c+k = A*k(-1)^alpha+(1-delta)*k(-1);
end;

// setting value for exogenous variable
// providing exact value for steady state of
// endogenous variables
steady_state_model;
k = ((1-beta*(1-delta))/(beta*alpha*A))^(1/(alpha-1));
c = A*k^alpha-delta*k;
end;

initval;
A = 1;
end;

// display the steady-state
steady;

// initial conditionls
histval;
k(0) = 0.5*((1-beta*(1-delta))/(beta*alpha))^(1/(alpha-1));
end;

// computing the simulation
perfect_foresight_setup(periods=200);
perfect_foresight_solver;
// ploting results
plot(simulation("k"))
plot(simulation("c"))
