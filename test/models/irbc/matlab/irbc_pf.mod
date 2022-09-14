@#for j in 1:N
  var lc_@{j} ll_@{j} lk_@{j} i_@{j} la_@{j} y_@{j} tn_@{j};
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
    y_@{j} = A*exp(alpha*lk_@{j}(-1) +(1 - alpha)*ll_@{j});
    t_@{j}*exp(lc_@{j})^(-1/gamma_@{j}) = lambda;
    t_@{j}*b_@{j}*exp(ll_@{j})^(1/eta_@{j}) = lambda*exp(la_@{j})*(1 - alpha)*y_@{j}/exp(ll_@{j});
    lambda*(1 +phi*(i_@{j}/exp(lk_@{j}(-1)) - delta))
      = beta*lambda(+1)*(1 + exp(la_@{j}(+1))*alpha*y_@{j}(+1)/exp(lk_@{j})
        + phi*(1 - delta + i_@{j}(+1)/exp(lk_@{j}) - (1/2)*(i_@{j}(+1)/exp(lk_@{j}) - delta))
             *(i_@{j}(+1)/exp(lk_@{j}) - delta));
    exp(lk_@{j}) = (1 - delta)*exp(lk_@{j}(-1)) + i_@{j};
    la_@{j} = rho*la_@{j}(-1) + sigma*(e + e_@{j});
  @#endfor
  tn_1 = exp(lc_1) + i_1 - delta*exp(lk_1(-1)) - y_1 + (phi/2)*exp(lk_1(-1))*(i_1/exp(lk_1(-1)) - delta)^2;
  @#for j in 2:N
  tn_@{j} =  tn_@{j-1} + exp(lc_@{j}) + i_@{j} - delta*exp(lk_@{j}(-1)) - y_@{j}
             + (phi/2)*exp(lk_@{j}(-1))*(i_@{j}/exp(lk_@{j}(-1)) - delta)^2;
  @#endfor
  tn_@{N} = 0;
end;


steady_state_model;
  @#for j in 1:N
    ll_@{j} = 0;
    lk_@{j} = 0;
    la_@{j} = 0;
    i_@{j} = delta;
    lc_@{j} = log(A);
    y_@{j} = A;
    tn_@{j} = 0;
  @#endfor
  lambda = 1;
end;

shocks;
  var e;
  periods 1;
  values 20;
  @#for j in 1:N
  var e_@{j};
  periods 1;
  values (2*@{j}/@{N});
  @#endfor
end;

perfect_foresight_setup(periods = 300);
perfect_foresight_solver;

