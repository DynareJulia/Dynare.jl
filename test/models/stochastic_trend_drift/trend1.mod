var y ygap ypot infl r;
varexo eg ep;

parameters g rho phi gamma theta;

rho = 0.8;
phi = 0.1;
gamma = 1.5;
theta = 0.1;
g = 2;

model;
  y = ygap + ypot;
  ypot = ypot(-1) + g + 0.2*ep;
  ygap = rho*ygap(-1) - phi*r + eg;
  infl = infl(+1) + theta*ygap;
  r = gamma*infl;
end;

shocks;
  var eg; stderr 1.0;
  var ep; stderr 1.0;
end;

deterministic_trends;
  y (g);
  ypot (g);
end;

histval;
  y(0) = 100;
  ypot(0) = 100;
  ygap(0) = 0;
end;

steady(nocheck);

stoch_simul(order=1, irf=0, nonstationary);

//CSV.write("data.csv", getfield(context.results.model_results[1].simulations[1].data, :data));

//histval_file(datafile=history);

//stoch_simul(periods=100, irf=0);
