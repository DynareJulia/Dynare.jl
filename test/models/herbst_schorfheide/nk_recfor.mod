var y R pi z g YGR INFL INT;
varexo eg ep ez;

parameters tau rho_R psi1 psi2 rho_g rho_z kappa  piA gammaQ rA s_ep s_eg s_ez; 

rA = 0.42;
piA = 3.3;
gammaQ = 0.52;
tau = 2.83;
kappa = 0.78;
psi1 = 1.8;
psi2 = 0.63;
rho_R = 0.77;
rho_g = 0.98;
rho_z = 0.88;
s_ep = 0.22;
s_eg = 0.72;
s_ez = 0.31;


model;
#  beta = 1/(1 + rA/400);
  y = y(+1) - (1/tau)*(R - pi(+1) - z(+1)) + g - g(+1);
  pi = beta*pi(+1) + kappa*(y - g);
  R = rho_R*R(-1) + (1 - rho_R)*(psi1*pi + psi2*(y - g)) + s_ep*ep/100;
  g = rho_g*g(-1) + s_eg*eg/100;
  z = rho_z*z(-1) + s_ez*ez/100;
  YGR = gammaQ + 100*(y - y(-1) + z);
  INFL = piA + 400*pi;
  INT = piA + rA + 4*gammaQ + 400*R;
end;

steady_state_model;
  y = 0;
  pi = 0;
  R = 0;
  g = 0;
  z = 0;
  YGR = gammaQ;
  INFL = piA;
  INT = piA + rA + 4*gammaQ;
end;

varobs YGR INFL INT;

shocks;
  var ep; stderr 1.0;
  var eg; stderr 1.0;
  var ez; stderr 1.0;
  // calibrated measurement errors
  var YGR;  stderr (0.20*0.579923);
  var INFL; stderr (0.20*1.470832);
  var INT;  stderr(0.20*2.237937);
end;

stoch_simul(order=1, irf=0);

recursive_forecasting!(periods=12, first_period=78, last_period=80, datafile="dsge1_data.csv");


