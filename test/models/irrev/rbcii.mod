var k, y, L, c, i, A, a, mu;
varexo epsilon;
parameters beta, theta, tau, alpha, psi, delta, rho, Astar, sigma;

beta    =  0.990;
theta   =  0.357;
tau     =  2.000;
alpha   =  0.450;
psi     = -2.500;
delta   =  0.020;
rho     =  0.998;
Astar   =  1.000;
sigma   =  0.100;

model;
 a = rho*a(-1) + sigma*epsilon;
 A = Astar*exp(a);
 (c^theta*(1-L)^(1-theta))^(1-tau)/c - mu = beta*((c(+1)^theta*(1-L(+1))^(1-theta))^(1-tau)/c(+1)*(alpha*(y(+1)/k)^(1-psi)+1-delta)-mu(+1)*(1-delta));
 ((1-theta)/theta)*(c/(1-L)) - (1-alpha)*(y/L)^(1-psi);
 y = A*(alpha*(k(-1)^psi)+(1-alpha)*(L^psi))^(1/psi);
 k = y-c+(1-delta)*k(-1);
 i = k-(1-delta)*k(-1);

[ mcp = 'i > 0' ]
 mu = 0;
end;

steady_state_model;
 a=0;
 mu=0;
 A=Astar;

 // Steady state ratios
 Output_per_unit_of_Capital=((1/beta-1+delta)/alpha)^(1/(1-psi));
 Consumption_per_unit_of_Capital=Output_per_unit_of_Capital-delta;
 Labour_per_unit_of_Capital=(((Output_per_unit_of_Capital/A)^psi-alpha)/(1-alpha))^(1/psi);
 Output_per_unit_of_Labour=Output_per_unit_of_Capital/Labour_per_unit_of_Capital;
 Consumption_per_unit_of_Labour=Consumption_per_unit_of_Capital/Labour_per_unit_of_Capital;

 L=1/(1+Consumption_per_unit_of_Labour/((1-alpha)*theta/(1-theta)*Output_per_unit_of_Labour^(1-psi)));
 c=Consumption_per_unit_of_Labour*L;
 k=L/Labour_per_unit_of_Capital;
 y=Output_per_unit_of_Capital*k;
 i=delta*k;
end;

shocks;
 var epsilon;
 periods 10;
 values -1;
end;

steady;

perfect_foresight_setup(periods=400);
perfect_foresight_solver(lmmcp, maxit=200, no_homotopy);
