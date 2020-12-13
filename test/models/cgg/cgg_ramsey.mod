var y pie r;
varexo e_y e_pie;

parameters delta sigma alpha kappa;

delta =  0.44;
kappa =  0.18;
alpha =  0.48;
sigma = -0.06;

model(linear);
y  = delta * y(-1)  + (1-delta)*y(+1)+sigma *(r - pie(+1)) + e_y; 
pie  =   alpha * pie(-1) + (1-alpha) * pie(+1) + kappa*y + e_pie;
end;

shocks;
var e_y;
stderr 0.63;
var e_pie;
stderr 0.4;
end;

planner_objective pie^2 + y^2;

ramsey_model(planner_discount=1.0);

stoch_simul(irf=0);

