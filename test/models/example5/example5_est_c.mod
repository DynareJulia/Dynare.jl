using CSV
//using Random
using TimeDataFrames

//Random.seed!(11)

var y, c, lk, a, lh, b;
varexo e, u;

parameters beta, rho, alpha, delta, theta, tau;

alpha = 0.36;
rho = 0.95;
tau = 0.025;
beta = 0.99;
delta = 0.025;
theta = 2.95;

model;
c*theta*exp(lh)=(1-alpha)*y;
exp(lk) = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
  *(exp(b(+1))*alpha*y(+1)+(1-delta)*exp(lk)));
y = exp(a)*(exp(lk(-1))^alpha)*(exp(lh)^(1-alpha));
exp(lk) = exp(b)*(y-c)+(1-delta)*exp(lk(-1));
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

steady_state_model;
  K_Y = beta*alpha /(1 - beta*(1 - delta));
  H_Y = K_Y^(-alpha/(1 - alpha));
  C_Y = 1 - delta*K_Y;
  y = (1 - alpha)/(theta*C_Y*H_Y);
  c = C_Y*y;
  lk = log(K_Y*y);
  lh = log(H_Y*y);
  a = 0;
  b = 0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
end;

stoch_simul(order=1, periods=100, irf=0);
CSV.write("fsdata_simul.csv", dataframe(simulation(1)))

//alpha.prior(shape=beta, mean=0.356, stdev=0.02);
//beta.prior(shape=beta, mean=0.993, stdev=0.002);
//rho.prior(shape=beta, mean=0.129, stdev=0.223);
//delta.prior(shape=beta, mean=0.01, stdev=0.005);
//theta.prior(shape=normal, mean=3.0, stdev=1.0);
//tau.prior(shape=beta, mean=0.03, stdev=0.01);
//std(e).prior(shape=inv_gamma, mean=0.009, stdev=100);
//std(u).prior(shape=inv_gamma, mean=0.009, stdev=100);

alpha.prior(shape=beta, mean=0.356, stdev=0.02);
beta.prior(shape=beta, mean=0.9, stdev=0.01);
rho.prior(shape=beta, mean=0.129, stdev=0.01);
delta.prior(shape=beta, mean=0.01, stdev=0.005);
theta.prior(shape=normal, mean=3.0, stdev=0.1);
tau.prior(shape=beta, mean=0.03, stdev=0.01);
std(e).prior(shape=inv_gamma, mean=0.009, stdev=1);
std(u).prior(shape=inv_gamma, mean=0.009, stdev=1);

varobs y, c;
//estimation(order=1,datafile=fsdata_simul);
//ml_estimation(context; datafile="fsdata_simul.csv");
//datafile = "fsdata_simul.csv";
hmc_estimation(context, datafile="fsdata_simul.csv");

