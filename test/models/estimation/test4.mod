// Inspired by Herbst and Schorfheide models
// Priors trighten in order to facilitate estimation
// Use prior!() syntax
// mode_compute! and rwmh_compute!
// estimated_params_init
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
# beta = 1/(1 + rA/400);
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

steady;

shocks;
  var ep; stderr 1.0;
  var eg; stderr 1.0;
  var ez; stderr 1.0;
  // calibrated measurement errors
  var YGR;  stderr (0.20*0.579923);
  var INFL; stderr (0.20*1.470832);
  var INT;  stderr(0.20*2.237937);
end;

estimated_params;
rA, gamma_pdf, 0.5, 0.2;
piA, gamma_pdf, 7, 2;
gammaQ, normal_pdf, 0.4, 0.2;
tau, gamma_pdf, 2, 0.5;
kappa, beta_pdf, 0.5, 0.2;
psi1, gamma_pdf, 1.5, 0.25;
psi2, gamma_pdf, 0.5, 0.25;
rho_R, beta_pdf, 0.5, 0.2;
rho_g, beta_pdf, 0.5, 0.2;
rho_z, beta_pdf, 0.5, 0.2;
s_eg, inv_gamma1_pdf, 1.2533141373155003, sqrt(0.4292036732051032);
s_ez, inv_gamma1_pdf, 0.6266570686577502, sqrt(0.1073009183012758);

stderr ep, inv_gamma1_pdf, 1.0, 0.1;
stderr INFL, inv_gamma1_pdf, 0.4, 0.1;

corr eg, ez, normal_pdf, 0.0, 0.01;
corr YGR, INT, normal_pdf, 0.0, 0.01;
end;

estimated_params_init;
  tau, 2.263895645223946;
  kappa, 0.9;
  psi1, 1.932357745633903;
  psi2, 0.465714640873466;
  rA, 0.333453026636814;
  piA, 3.427852381637088;
  gammaQ, 0.624767199308368;
  rho_R, 0.764476492352786;
  rho_g, 0.990359718181570;
  rho_z, 0.917302546854851;
  s_ez, 0.192160214260411;
  s_eg, 0.632177304244290;
end;

varobs YGR INFL INT;

estimation(datafile='dsge1_data.csv', mh_replic=100000, mh_jscale=0.15, mh_nblocks=2);
