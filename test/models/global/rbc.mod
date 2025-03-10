var k z y c rk;

varexo e;

parameters alpha beta theta rho sig
           ss_l ss_k ss_c ss_y ss_rk;

alpha   = 1/3;
beta    = 0.99;
rho     = 0.95;
sig     = 0.01;
ss_l    = 1/3;
ss_k    = ss_l*(beta*alpha)^(1/(1-alpha));
ss_y    = ss_k^alpha*ss_l^(1-alpha);
ss_c    = ss_y-ss_k;
ss_rk   = alpha*ss_y/ss_k;
theta   = (1-alpha)*(1-ss_l)*ss_y/(ss_l*ss_c);

model; 
  z   = rho*z(-1)+sig*e;
  1/c = beta*rk(+1)/c(+1);
  y   = exp(z)*k(-1)^alpha*ss_l^(1-alpha);
  rk  = alpha*y/k(-1);
  y   = c+k;
end;

shocks;
  var e; stderr 1;
end;

initval;
  k = ss_k;
  z = 0;
  c = ss_c;
  rk = ss_rk;
  e = 0;
end;

steady;

limits!("k", min = 0.1*ss_k, max = 1.9*ss_k);
limits!("z", min = -2*sig/sqrt(1-rho^2), max = 2*sig/sqrt(1-rho^2));

