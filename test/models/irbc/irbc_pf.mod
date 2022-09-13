//@#define N=2
@#for j in 1:N
  var c_@{j} l_@{j} k_@{j} i_@{j} a_@{j} y_@{j} tn_@{j};
  varexo e_@{j};
  parameters b_@{j} gamma_@{j} eta_@{j} t_@{j};

  gamma_@{j} = 0.25;
  eta_@{j}   = 0.1;  

@#endfor

var lambda;
varexo e;
parameters alpha beta delta phi rho A sigma;

alpha = 0.36;
beta  = 0.99;  
delta = 0.025;
phi   = 0.5;
rho   = 0.95;
sigma = 0.01;
A = (1 - beta)/(alpha*beta);
@#for j in 1:N
  t_@{j} = 1/(A^(-1/gamma_@{j}));
  b_@{j} = (1 - alpha)*A^(1-1/gamma_@{j});
@#endfor

model;
  @#for j in 1:N
    y_@{j} = a_@{j}*A*k_@{j}(-1)^alpha*l_@{j}^(1 - alpha);
    t_@{j}*c_@{j}^(-1/gamma_@{j}) = lambda;
    t_@{j}*b_@{j}*l_@{j}^(1/eta_@{j}) = lambda*a_@{j}*(1 - alpha)*y_@{j}/l_@{j};
    lambda*(1 +phi*(i_@{j}/k_@{j}(-1) - delta))
      = beta*lambda(+1)*(1 + a_@{j}(+1)*alpha*y_@{j}(+1)/k_@{j}
        + phi*(1 - delta + i_@{j}(+1)/k_@{j} - (1/2)*(i_@{j}(+1)/k_@{j} - delta))
             *(i_@{j}(+1)/k_@{j} - delta));
    k_@{j} = (1 - delta)*k_@{j}(-1) + i_@{j};
    log(a_@{j}) = rho*log(a_@{j}(-1)) + sigma*(e + e_@{j});
  @#endfor
  tn_1 = c_1 + i_1 - delta*k_1(-1) - y_1 + (phi/2)*k_1(-1)*(i_1/k_1(-1) - delta)^2;
  @#for j in 2:N
  tn_@{j} =  tn_@{j-1} + c_@{j} + i_@{j} - delta*k_@{j}(-1) - y_@{j}
             + (phi/2)*k_@{j}(-1)*(i_@{j}/k_@{j}(-1) - delta)^2;
  @#endfor
  tn_@{N} = 0;
end;


steady_state_model;
  @#for j in 1:N
    l_@{j} = 1;
    k_@{j} = 1;
    a_@{j} = 1;
    i_@{j} = delta;
    c_@{j} = A;
    y_@{j} = A;
    tn_@{j} = 0;
  @#endfor
  lambda = 1;
end;

shocks;
  var e; stderr 0.01;
  @#for j in 1:N
    var e_@{j}; stderr 0.01;
  @#endfor
end;

check;

stoch_simul(order=1, irf=0);

