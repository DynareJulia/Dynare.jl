///////
// IRBC model from
// Brumm, Krause, Schaab, Scheidegger "Sparse grids for dynamic economic models", 2021.
// https://github.com/SparseGridsForDynamicEcon/SparseGrids_in_econ_handbook/blob/master/doc/sparse_grids_in_econ.pdf
//
// Utility function:
//   u_j(c_j_t) = c_j_t^(1-gamma_j)/(1 - gamma_j)
// Production function:
//   a_j_t * A_t * k_j_{t-1}^kappa 
// Capital adjustment cost:
//   0.5 * phi * k_j_{t-1} * (k_j_t/k_j_{t-1} - 1)^2
//
///////

@#define N=2
@#for j in 1:N
var lk_@{j} a_@{j} i_@{j};
  varexo e_@{j};
  parameters gamma_@{j} eta_@{j} t_@{j};

  a_eis = 0.25;
  b_eis = 1;
  gamma_@{j} = a_eis + (@{j} - 1)*(b_eis - a_eis)/(@{N}-1);
  eta_@{j} = 0.1;  

@#endfor

var llambda;
varexo e;
parameters kappa beta delta phi rho A sigE;

// zeta
kappa = 0.36;
// betta
beta  = 0.99;  
delta = 0.01;
// kappa
phi   = 0.5;
// rhoz
rho   = 0.95;
// sigE
sigE = 0.01;
// A_tfp
A = (1 - beta*(1 - delta))/(kappa*beta);
@#for j in 1:N
  // pareto
//  t_@{j} = (A - delta)^(1/gamma_@{j});
  t_@{j} = A^(1/gamma_@{j});
@#endfor

model;
  @#for j in 1:N
    i_@{j} = exp(lk_@{j}) - (1 - delta)*exp(lk_@{j}(-1));
    [mcp = 'i_@{j} > 0']
    exp(llambda)*(1 + phi*(exp(lk_@{j} - lk_@{j}(-1)) - 1))
      = beta*exp(llambda(+1))*(exp(a_@{j}(+1))*kappa*A*exp((kappa-1)*lk_@{j})
        + 1 - delta + (phi/2)*(exp(lk_@{j}(+1) - lk_@{j}) - 1)*(exp(lk_@{j}(+1) - lk_@{j}) + 1));
   [preamble]
      a_@{j} = rho*a_@{j}(-1) + sigE*(e + e_@{j});
  @#endfor
    exp(a_1)*A*exp(kappa*lk_1(-1))
  @#for j in 2:N
    + exp(a_@{j})*A*exp(kappa*lk_@{j}(-1))
  @#endfor
  =
  (exp(llambda)/t_1)^(-gamma_1) + exp(lk_1) - (1 - delta)*exp(lk_1(-1)) + (phi/2)*exp(lk_1(-1))*(exp(lk_1 - lk_1(-1)) - 1)^2
  @#for j in 2:N
  + (exp(llambda)/t_@{j})^(-gamma_@{j}) + exp(lk_@{j}) - (1 - delta)*exp(lk_@{j}(-1))
             + (phi/2)*exp(lk_@{j}(-1))*(exp(lk_@{j} - lk_@{j}(-1)) - 1)^2
  @#endfor
  ;
end;


//steady_state_model;
initval;
@#for j in 1:N
    lk_@{j} = 0;
    a_@{j} = 0;
    i_@{j} = delta;
  @#endfor
  llambda = 0;
end;

steady;

shocks;
  var e; stderr 1;
  @#for j in 1:N
    var e_@{j}; stderr 1;
  @#endfor
end;

//check;

//stoch_simul(order=1, irf=0);

@#for j in 1:N
  limits!("lk_@{j}", min = log(0.8), max = log(1.2));
  limits!("a_@{j}", min = -0.8*sigE/(1 - rho), max = 0.8*sigE/(1 - rho));
@#endfor

//shocks;
//  var e_1;
//  periods 1:10;
//  values -4;
//end;

(grid, state_variables, policy_variables) = sparsegridapproximation(scaleCorrExclude=["lambda"], mcp = true);

Y = simulate!(context, grid, 1000000, policy_variables, state_variables, 1);
M = mean(Y, dims=3);
display(Plots.plot(M[:, [1,4]], label=["lk_1" "lk_2"], title="10-period shock on e_1"));
