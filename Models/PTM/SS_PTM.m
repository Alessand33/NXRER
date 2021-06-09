%============================================================
%  Deterministic Steady State Computation
%============================================================


function F = SS_PTM(x,xpar);

bet    = xpar(1);
sig    = xpar(2); 
th     = xpar(3);
rho    = xpar(4);
sdeta  = xpar(5);
PsiT   = xpar(6); 
Psi0   = xpar(7);  
Psi1   = xpar(8); 
PsiX   = xpar(9); 
n1     = xpar(10); 
n0     = xpar(11); 
N      = xpar(12); 
eta0   = xpar(13); 
eta1   = xpar(14); 
Tr     = xpar(15); 
L      = xpar(16); 
sdb    = xpar(17);


C    = exp(x(1)); 
tau0 = exp(x(2)); 
tau1 = exp(x(3));
a2   = exp(x(4)); 
gam  = exp(x(5));
dEV  = exp(x(6));
Ph   = exp(x(7));
Pf   = exp(x(8));


W = (1-gam)/gam*C/(1-L);
Lp  = L - (1-N)*n0*tau0 - N*(1-n1)*tau1;
EXY =  a2*Pf^(1-rho);

F(1) = 1/th*a2*(th*W/(th-1))^(1-th)*exp((th-1)*eta0)*Pf^(th-rho)*C + bet*(dEV) - W*tau0;

F(2) = 1/th*a2*(th*W/(th-1))^(1-th)*exp((th-1)*eta1)*Pf^(th-rho)*C + bet*(dEV) - W*tau1;
    
F(3) = (th-1)/th*C/W - Lp;

F(4) = (th*W/(th-1))^(1-th)*PsiT - Ph^(1-th);

F(5) = (th*W/(th-1))^(1-th)*PsiX - Pf^(1-th);

F(6) = 1/th*a2*(th*W/(th-1))^(1-th)*(Psi1-Psi0)*Pf^(th-rho)*C - ((1-n1)*W*tau1-n0*W*tau0) + bet*((1-n1)-n0)*(dEV) - dEV;
 
F(7) = Tr - EXY;

F(8) = Ph^(1-rho) + a2*Pf^(1-rho) -1;


F=F';

