// test case for reporting
using CSV

var y, c, k, a, h, b;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*0.5*(c(+1)+c(+2))))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

steady_state_model;
  K_Y = beta*alpha/(1 - beta*(1 - delta));
  H_Y = K_Y^(-alpha/(1 - alpha));
  C_Y = 1 - delta*K_Y;
  y = (theta*C_Y*H_Y^(1 + psi)/(1 - alpha))^(-1/(1 + psi));
  c = C_Y*y;
  k = K_Y*y;
  h = H_Y*y;
  a = 0;
  b = 0;
end;

steady;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = 0.1*0.009*0.009;
end;

stoch_simul(order=1, periods=100);

CSV.write("data.csv", simulation(1))

varobs y;
calib_smoother(datafile='data.csv');

// Reporting
report = Report("Report example", subtitle="Example3report.mod")

page1 = Page()
page2 = Page()
page3 = Page()

res = context.results.model_results[1]
res1 = Matrix{Any}(permutedims(res.linearrationalexpectations.g1[1:6,1:6], [2, 1]))
table = Table(res1, "Reduced form", ["y" "c" "k" "a" "h" "b"], LatexCell.(["\$ \\phi(c) \$", "\$ \\phi(k) \$", "\$ \\phi(a)\$", "\$ \\phi(b) \$", "e", "u"]), "Note: \$ \\phi(x) = x_{t-1} - steady\\_state(x)\$")
 
add_model!(page1, lastline = 52)

add_table!(page2, table)

plot(res.smoother[:y], filename="graph.png")

add_graph!(page3, Graph("graph.png"))

add_page!(report, page1)
add_page!(report, page2)
add_page!(report, page3)

print(report)
