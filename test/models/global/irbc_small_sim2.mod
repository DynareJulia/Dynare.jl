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
var k_@{j} a_@{j};
  varexo e_@{j};
  parameters gamma_@{j} eta_@{j} t_@{j};

  a_eis = 0.25;
  b_eis = 1;
  gamma_@{j} = a_eis + (@{j} - 1)*(b_eis - a_eis)/(@{N}-1);
  eta_@{j} = 0.1;  

@#endfor

var lambda;
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
  t_@{j} = A^(1/gamma_@{j});
@#endfor

model;
  @#for j in 1:N
    lambda*(1 + phi*(k_@{j}/k_@{j}(-1) - 1))
      = beta*lambda(+1)*(exp(a_@{j}(+1))*kappa*A*k_@{j}^(kappa - 1)
        + 1 - delta + (phi/2)*(k_@{j}(+1)/k_@{j} - 1)*(k_@{j}(+1)/k_@{j} + 1));
   [preamble]
      a_@{j} = rho*a_@{j}(-1) + sigE*(e + e_@{j});
  @#endfor
    exp(a_1)*A*k_1(-1)^kappa
  @#for j in 2:N
    + exp(a_@{j})*A*k_@{j}(-1)^kappa
  @#endfor
  =
  (lambda/t_1)^(-gamma_1) + k_1 - (1 - delta)*k_1(-1) + (phi/2)*k_1(-1)*(k_1/k_1(-1) - 1)^2
  @#for j in 2:N
  + (lambda/t_@{j})^(-gamma_@{j}) + k_@{j} - (1 - delta)*k_@{j}(-1)
             + (phi/2)*k_@{j}(-1)*(k_@{j}/k_@{j}(-1) - 1)^2
  @#endfor
  ;
end;


//steady_state_model;
initval;
@#for j in 1:N
    k_@{j} = 1;
    a_@{j} = 0;
  @#endfor
  lambda = 1;
end;

steady;

shocks;
  var e; stderr 1;
  @#for j in 1:N
    var e_@{j}; stderr 1;
  @#endfor
end;

@#for j in 1:N
  limits!("k_@{j}", min = 0.8, max = 1.2);
  limits!("a_@{j}", min = -0.8*sigE/(1 - rho), max = 0.8*sigE/(1 - rho));
@#endfor

(grid, sgws) = sparsegridapproximation(scaleCorrExclude=["lambda"]);
Y = simulate!();
errors = simulation_approximation_error();

// Y = simulate!(context, grid, 100, policy_variables, state_variables, 10000);
// M = mean(Y, dims=3);
// display(Plots.plot(M[:, [1,3]], label=["k_1" "k_2"], title="10-period shock on e_1"));

