//MODEL:
//test on Dynare to find the lagrangean multipliers.
//We consider a standard NK model. We use the FOCS of the competitive economy and we aim at calculating the Ramsey optimal problem.

//------------------------------------------------------------------------------------------------------------------------
//1. Variable declaration
//------------------------------------------------------------------------------------------------------------------------

var pai, c, n, r, a;

//4 variables + 1 shock

varexo u;




//------------------------------------------------------------------------------------------------------------------------
// 2. Parameter declaration and calibration
//-------------------------------------------------------------------------------------------------------------------------

parameters beta, rho, epsilon, omega, phi, gamma;

beta=0.99;
gamma=3;       //Frish elasticity
omega=17;        //price stickyness
epsilon=8;      //elasticity for each variety of consumption
phi=1;           //coefficient associated to labor effort disutility

rho=0.95;  //coefficient associated to productivity shock


//-----------------------------------------------------------------------------------------------------------------------
// 3. The model
//-----------------------------------------------------------------------------------------------------------------------


model;


a=rho*(a(-1))+u;

1/c=beta*(1/(c(+1)))*(r/(pai(+1)));               //euler


omega*pai*(pai-1)=beta*omega*(c/(c(+1)))*(pai(+1))*(pai(+1)-1)+epsilon*exp(a)*n*(c/exp(a)*phi*n^gamma-(epsilon-1)/epsilon);  //NK pc
//pai*(pai-1)/c = beta*pai(+1)*(pai(+1)-1)/c(+1)+epsilon*phi*n^(gamma+1)/omega-exp(a)*n*(epsilon-1)/(omega*c);  //NK pc

(exp(a))*n=c+(omega/2)*((pai-1)^2);

end;

//--------------------------------------------------------------------------------------------------------------------------
// 4. Steady state
//---------------------------------------------------------------------------------------------------------------------------

initval;
pai=1;
r=1;
c=11;
n=1;
a=0;
end;



//---------------------------------------------------------------------------------------------------------------------------
// 5. shocks
//---------------------------------------------------------------------------------------------------------------------------

shocks;
var u; stderr 0.008;

end;

planner_objective(ln(c)-phi*((n^(1+gamma))/(1+gamma)));

ramsey_model(planner_discount=0.99);

steady;
