%============================================================
%  Deterministic Steady State Computation
%============================================================


function F = SteadyS(x);

load xpar;
load xss0;

bet    = xpar(1);
sig    = xpar(2); 
th     = xpar(3);
els    = xpar(4);
rho    = xpar(5);
alp    = xpar(6);
del    = xpar(7);
sdeta  = xpar(8);
rhoz   = xpar(9);
sdeps  = xpar(10);
coreps = xpar(11);
co     = xpar(12);
mu     = xpar(13);
v      = xpar(14);
m      = xpar(15);
aCeqaI = xpar(16);

Tr     = xss0(1); 
L      = xss0(2); 

C    = exp(x(1)); 
Ph   = exp(x(2)); 
Pf   = exp(x(3)); 
K1   = exp(x(4)); 
a2   = exp(x(5)); 
gam  = exp(x(6));
a2I  = exp(x(7));

W = (1-gam)/gam*C/(1-L);
a1 = 1;
K  = K1;
Iv  = del*K;
Lp   = L;
PI  = ( Ph^(rho/(rho-1)) + a2I^(1/(1-rho))*Pf^(rho/(rho-1)) )^((rho-1)/rho);
Dh = C + PI^(1/(1-rho))*Iv;
Df = a2^(1/(1-rho))*C + a2I^(1/(1-rho))*PI^(1/(1-rho))*Iv;
L1 = co^(v/(1-th))*W^(v/(th-1))*K1^(1-v)*( ( Ph^mu*Dh + Pf^mu*Df )^v );
EXN = Pf^(rho/(rho-1))*Df;
YN = 1/co*W*Lp;
      
F(1) = co^(v/(1-th)-1)*W^(1+v/(th-1))*( Ph^mu*Dh + Pf^mu*Df )^(v-1)*( K1^(1-v) )...
        - Pf^(th/(th-1));

F(2) = Pf^(th/(th-1))...
        - Ph^(th/(th-1));
    
F(3) = Ph^(rho/(rho-1)) + a2^(1/(1-rho))*Pf^(rho/(rho-1)) - 1;

F(4) = bet*( alp/(1-alp)*W*L1/(PI*K1) + 1-del) - 1;

F(5) = EXN/YN - Tr;

F(6) = L1 - Lp;

%F(10) = a2I - 1;
F(7) = a2I - (aCeqaI<0.5) - (aCeqaI>0.5)*a2;




F=F';

