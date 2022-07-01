// one auxiliary variable
// Generalized Schur algorithm

var y, c, k, a, h, b;
varexo e, u;

parameters beta, rho, alpha, delta, theta, tau;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
theta = 2.95;

model;
c*theta*h^(1)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*0.5*(c(+1)+c(+2))))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

steady_state_model;
  K_Y = beta*alpha /(1 - beta*(1 - delta));
  H_Y = K_Y^(-alpha/(1 - alpha));
  C_Y = 1 - delta*K_Y;
  y = (theta*C_Y*H_Y^(1 + psi)/(1 - alpha))^(-1/(1 + psi));
  c = C_Y*y;
  k = K_Y*y;
  h = H_Y*y;
  a = 0;
  b = 0;
end;
  
shocks;
var e; stderr 0.009;
var u; stderr 0.009;
end;

stoch_simul(order=1, periods=100, irf=0);



