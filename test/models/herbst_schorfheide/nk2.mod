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

shocks;
  var ep; stderr 1.0;
  var eg; stderr 1.0;
  var ez; stderr 1.0;
end;

estimated_params;
rA, gamma_pdf, 0.5, 0.5;
piA, gamma_pdf, 7, 2;
gammaQ, normal_pdf, 0.4, 0.2;
tau, gamma_pdf, 2, 0.5;
//kappa, uniform_pdf,,,0,1;
kappa, beta_pdf, 0.5, 0.1;
psi1, gamma_pdf, 1.5, 0.25;
psi2, gamma_pdf, 0.5, 0.25;
rho_R, beta_pdf, 0.5, 0.1;
rho_g, beta_pdf, 0.5, 0.15;
rho_z, beta_pdf, 0.5, 0.1;
// s=0.4, nu=4
s_ep, inv_gamma2_pdf, 0.32, inf;
// s=1, nu=4
s_eg, inv_gamma2_pdf, 1.0, inf;
// s=0.5, nu=4
s_ez, inv_gamma2_pdf, 0.5, inf;
end;

//stoch_simul(order=1, irf=0);

estimated_params_init;
rA, 0.42;
piA, 3.3;
gammaQ, 0.52;
tau, 2.83;
kappa, 0.78;
psi1, 1.8;
psi2, 0.63;
rho_R, 0.77;
rho_g, 0.98;
rho_z, 0.88;
s_ep, 0.22;
s_eg, 0.71;
s_ez, 0.31;
end;

varobs YGR INFL INT;

estimation(datafile='dsge1_data.csv', mh_replic= 0, mode_check);
