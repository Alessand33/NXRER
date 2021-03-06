function [residual, g1, g2, g3] = dynare_NoK02_decomp_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(74, 1);
T10 = (1-params(18))/params(18);
T32 = exp(y(62));
T35 = 1/T32*params(17);
T39 = exp(y(58))^(1-T32);
T43 = T32-1;
T44 = T32*exp(y(27))/T43;
T45 = T44^(1-T32);
T46 = T35*T39*T45;
T48 = params(13)^2;
T52 = exp(T48*T43^2/2);
T53 = T46*T52;
T57 = 1-normcdf(y(33),T43*T48,params(13));
T70 = exp(y(14)+T32*y(19)+T43*y(59)+(T32-params(4))*y(39));
T71 = T53*T57*T70;
T90 = exp(params(2)*(y(13)-y(85))+(1-params(18))*(1-params(2))*(y(27)-y(88)));
T91 = exp(y(25))*T90;
T99 = exp(y(90))+exp(y(53))*(exp(y(92))-exp(y(90)));
T108 = params(17)*1/exp(y(61));
T112 = exp(y(57))^(1-exp(y(61)));
T116 = exp(y(61))-1;
T117 = exp(y(61))*exp(y(28))/T116;
T118 = T117^(1-exp(y(61)));
T119 = T108*T112*T118;
T123 = exp(T48*T116^2/2);
T124 = T119*T123;
T128 = 1-normcdf(y(34),T48*T116,params(13));
T140 = exp(y(13)+y(19)*(-exp(y(61)))+T116*y(60)+(exp(y(61))-params(4))*y(40));
T141 = T124*T128*T140;
T156 = exp(params(2)*(y(14)-y(86))+(1-params(18))*(1-params(2))*(y(28)-y(89)));
T165 = exp(y(91))+exp(y(54))*(exp(y(93))-exp(y(91)));
T166 = exp(y(26))*T156*T165;
T173 = 1-normcdf(y(35),T43*T48,params(13));
T184 = exp(y(90))+(exp(y(92))-exp(y(90)))*(1-exp(y(55)));
T192 = 1-normcdf(y(36),T48*T116,params(13));
T202 = exp(y(91))+(exp(y(93))-exp(y(91)))*(1-exp(y(56)));
T212 = exp(y(14)+(T32-params(4))*y(39)+T32*y(19)+T43*(y(33)+y(59)));
T223 = exp(y(13)+(exp(y(61))-params(4))*y(40)+y(19)*(-exp(y(61)))+T116*(y(34)+y(60)));
T234 = exp(y(14)+(T32-params(4))*y(39)+T32*y(19)+T43*(y(59)+y(35)));
T244 = exp(y(13)+(exp(y(61))-params(4))*y(40)+y(19)*(-exp(y(61)))+T116*(y(60)+y(36)));
T258 = exp(y(1))*(1-normcdf(y(35),0,params(13)));
T271 = exp(y(2))*(1-normcdf(y(36),0,params(13)));
T291 = exp(y(27))*params(3)/(params(3)-1);
T292 = T291^(-1);
T298 = T292*exp(y(13)+(1-params(4))*y(37));
T300 = params(17)*T44^(-1);
T304 = exp(y(14)+y(19)+y(39)*(1-params(4)));
T305 = T300*T304;
T309 = exp(y(28))*params(3)/(params(3)-1);
T310 = T309^(-1);
T315 = T310*exp(y(14)+(1-params(4))*y(38));
T317 = params(17)*T117^(-1);
T323 = T317*exp(y(13)+(-y(19))+y(40)*(1-params(4)));
T332 = params(3)*exp(y(27)-y(59))/(params(3)-1);
T342 = params(3)*exp(y(28)-y(60))/(params(3)-1);
T350 = T32*exp(y(27)-y(59))/T43;
T351 = T350^(1-T32);
T355 = T39*T351*exp(y(41));
T361 = exp(y(61))*exp(y(28)-y(60))/T116;
T362 = T361^(1-exp(y(61)));
T366 = T112*T362*exp(y(42));
T427 = 1+params(7)*y(44)*exp(y(43)-y(45));
T448 = exp((1-params(18))*(1-params(2))*(y(28)-y(89))+params(2)*(y(14)-y(86))+y(19)-y(87));
T449 = exp(y(26))*T448;
T455 = 1-params(7)*y(44)*exp(y(43)-y(46)-y(19));
T498 = exp(y(19)+y(39)-y(58)+y(23)*1/T43);
T511 = exp(y(40)-y(57)+y(24)*1/T116);
lhs =T10*exp(y(13)-y(27));
rhs =1-exp(y(15));
residual(1)= lhs-rhs;
lhs =T10*exp(y(14)-y(28));
rhs =1-exp(y(16));
residual(2)= lhs-rhs;
lhs =exp(y(29));
rhs =T71-exp(y(27)+y(53))*params(15)+T91*T99;
residual(3)= lhs-rhs;
lhs =exp(y(30));
rhs =T141-params(15)*exp(y(28)+y(54))+T166;
residual(4)= lhs-rhs;
lhs =exp(y(31));
rhs =T70*T53*T173-exp(y(27))*(1-exp(y(55)))*params(16)+T91*T184;
residual(5)= lhs-rhs;
lhs =exp(y(32));
rhs =T140*T124*T192-params(16)*exp(y(28))*(1-exp(y(56)))+exp(y(26))*T156*T202;
residual(6)= lhs-rhs;
lhs =exp(y(27))*params(15);
rhs =T46*T212+T91*(exp(y(92))-exp(y(90)));
residual(7)= lhs-rhs;
lhs =params(15)*exp(y(28));
rhs =T119*T223+exp(y(26))*T156*(exp(y(93))-exp(y(91)));
residual(8)= lhs-rhs;
lhs =exp(y(27))*params(16);
rhs =T91*(exp(y(92))-exp(y(90)))+T46*T234;
residual(9)= lhs-rhs;
lhs =exp(y(28))*params(16);
rhs =exp(y(26))*T156*(exp(y(93))-exp(y(91)))+T119*T244;
residual(10)= lhs-rhs;
lhs =exp(y(23));
rhs =(1-normcdf(y(33),0,params(13)))*(1-exp(y(1)))+T258;
residual(11)= lhs-rhs;
lhs =exp(y(24));
rhs =(1-normcdf(y(34),0,params(13)))*(1-exp(y(2)))+T271;
residual(12)= lhs-rhs;
lhs =exp(y(15));
rhs =exp(y(17))+params(15)*(1-normcdf(y(33),0,params(13)))*(1-exp(y(1)))+params(16)*T258;
residual(13)= lhs-rhs;
lhs =exp(y(16));
rhs =exp(y(18))+params(15)*(1-normcdf(y(34),0,params(13)))*(1-exp(y(2)))+params(16)*T271;
residual(14)= lhs-rhs;
lhs =exp(y(17));
rhs =T298+T305;
residual(15)= lhs-rhs;
lhs =exp(y(18));
rhs =T315+T323;
residual(16)= lhs-rhs;
lhs =exp(y(37)*(1-params(3)));
rhs =T332^(1-params(3))*params(14);
residual(17)= lhs-rhs;
lhs =exp(y(38)*(1-params(3)));
rhs =params(14)*T342^(1-params(3));
residual(18)= lhs-rhs;
lhs =exp((1-T32)*(y(19)+y(39)));
rhs =T355;
residual(19)= lhs-rhs;
lhs =exp((1-exp(y(61)))*(y(40)-y(19)));
rhs =T366;
residual(20)= lhs-rhs;
lhs =exp(y(41));
rhs =T173*T52*exp(y(1))+T57*T52*(1-exp(y(1)));
residual(21)= lhs-rhs;
lhs =exp(y(42));
rhs =T192*T123*exp(y(2))+T128*T123*(1-exp(y(2)));
residual(22)= lhs-rhs;
lhs =1;
rhs =exp((1-params(4))*y(37))+params(17)*exp(y(40)*(1-params(4)));
residual(23)= lhs-rhs;
lhs =1;
rhs =exp((1-params(4))*y(38))+params(17)*exp(y(39)*(1-params(4)));
residual(24)= lhs-rhs;
lhs =y(25);
rhs =(1-params(5))*log(params(1))+params(6)*log(params(19))+params(5)*y(3)-y(13)*params(6)-params(35)*x(it_, 5)/2;
residual(25)= lhs-rhs;
lhs =y(26);
rhs =params(35)*x(it_, 5)/2+(1-params(5))*log(params(1))+params(6)*log(params(19))+params(5)*y(4)-y(14)*params(6);
residual(26)= lhs-rhs;
lhs =exp(y(43))*T427;
rhs =T91;
residual(27)= lhs-rhs;
lhs =exp(y(13))+exp(y(43))*y(44);
rhs =exp(y(27)+y(17))+y(5)+exp(y(13)+(1-params(4))*y(37))*1/params(3)+T35*T304;
residual(28)= lhs-rhs;
lhs =T91/T427;
rhs =T449/T455;
residual(29)= lhs-rhs;
lhs =y(22);
rhs =(exp(y(47))-exp(y(48)))/exp(y(45));
residual(30)= lhs-rhs;
lhs =exp(y(45));
rhs =exp(y(13)+(1-params(4))*y(37))+params(17)*T304;
residual(31)= lhs-rhs;
lhs =exp(y(46));
rhs =exp(y(14)+(1-params(4))*y(38))+params(17)*exp(y(13)+(-y(19))+y(40)*(1-params(4)));
residual(32)= lhs-rhs;
lhs =exp(y(11));
rhs =exp(y(45)-y(37));
residual(33)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(46)-y(38));
residual(34)= lhs-rhs;
lhs =y(47);
rhs =y(14)+y(39)*(1-params(4))+y(19)+log(params(17));
residual(35)= lhs-rhs;
lhs =y(48);
rhs =y(13)+y(40)*(1-params(4))+log(params(17));
residual(36)= lhs-rhs;
lhs =exp(y(49));
rhs =T498+params(39)*x(it_, 6)/2;
residual(37)= lhs-rhs;
lhs =exp(y(50));
rhs =T511-params(39)*x(it_, 6)/2;
residual(38)= lhs-rhs;
lhs =exp(y(51));
rhs =exp(y(47)-y(49));
residual(39)= lhs-rhs;
lhs =exp(y(52));
rhs =exp(y(48)-y(50));
residual(40)= lhs-rhs;
lhs =1-normcdf(y(33),0,params(13));
rhs =exp(y(53));
residual(41)= lhs-rhs;
lhs =1-normcdf(y(34),0,params(13));
rhs =exp(y(54));
residual(42)= lhs-rhs;
lhs =normcdf(y(35),0,params(13));
rhs =exp(y(55));
residual(43)= lhs-rhs;
lhs =normcdf(y(36),0,params(13));
rhs =exp(y(56));
residual(44)= lhs-rhs;
lhs =y(57);
rhs =y(69)+0.5*y(70)+y(80);
residual(45)= lhs-rhs;
lhs =y(58);
rhs =y(69)-0.5*y(70);
residual(46)= lhs-rhs;
lhs =y(75);
rhs =params(8)*y(8)+params(9)*x(it_, 1);
residual(47)= lhs-rhs;
lhs =y(76);
rhs =params(23)*y(9)+params(22)*x(it_, 2);
residual(48)= lhs-rhs;
lhs =y(59);
rhs =y(75)+0.5*y(76);
residual(49)= lhs-rhs;
lhs =y(60);
rhs =y(75)-0.5*y(76);
residual(50)= lhs-rhs;
lhs =exp(y(61));
rhs =params(3)*exp(y(19)*params(12));
residual(51)= lhs-rhs;
lhs =T32;
rhs =params(3)*exp(y(19)*(-params(12)));
residual(52)= lhs-rhs;
lhs =y(20);
rhs =y(50)-y(49);
residual(53)= lhs-rhs;
lhs =y(21);
rhs =y(51)-y(52);
residual(54)= lhs-rhs;
lhs =y(63);
rhs =y(50)-y(49);
residual(55)= lhs-rhs;
lhs =y(64);
rhs =(exp(y(47))+exp(y(48)))/exp(y(45));
residual(56)= lhs-rhs;
lhs =y(65);
rhs =(exp(y(51))-exp(y(52)))/exp(y(11));
residual(57)= lhs-rhs;
lhs =y(66);
rhs =y(21)-(y(14)-y(13))-params(4)*y(72);
residual(58)= lhs-rhs;
lhs =y(67);
rhs =y(21)-(y(12)-y(11))-params(4)*y(72);
residual(59)= lhs-rhs;
lhs =y(68);
rhs =y(52)-y(11);
residual(60)= lhs-rhs;
lhs =y(69);
rhs =params(10)*y(6)+params(11)*x(it_, 3)/params(4);
residual(61)= lhs-rhs;
lhs =y(70);
rhs =params(20)*y(7)+params(21)*x(it_, 4)/params(4);
residual(62)= lhs-rhs;
lhs =y(71);
rhs =y(12)-y(11);
residual(63)= lhs-rhs;
lhs =y(72);
rhs =y(19)+y(20);
residual(64)= lhs-rhs;
lhs =y(73);
rhs =y(11)-params(24);
residual(65)= lhs-rhs;
lhs =y(74);
rhs =y(14)-y(13);
residual(66)= lhs-rhs;
lhs =y(77);
rhs =y(44)/exp(y(45));
residual(67)= lhs-rhs;
lhs =y(78);
rhs =y(47)-y(48);
residual(68)= lhs-rhs;
lhs =y(79);
rhs =y(72)+params(36);
residual(69)= lhs-rhs;
lhs =y(82);
rhs =y(74)+params(37);
residual(70)= lhs-rhs;
lhs =y(80);
rhs =params(38)*y(10);
residual(71)= lhs-rhs;
lhs =exp(y(81));
rhs =(exp(y(51))+exp(y(52)))/exp(y(11));
residual(72)= lhs-rhs;
lhs =y(83);
rhs =y(71)+params(37);
residual(73)= lhs-rhs;
lhs =y(84);
rhs =y(63);
residual(74)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(74, 99);

  %
  % Jacobian matrix
  %

T825 = T35*T39*T44*getPowerDeriv(T44,1-T32,1);
T826 = T52*T825;
T851 = getPowerDeriv(T44,(-1),1);
T857 = getPowerDeriv(T332,1-params(3),1);
T861 = getPowerDeriv(T350,1-T32,1);
T883 = T108*T112*T117*getPowerDeriv(T117,1-exp(y(61)),1);
T884 = T123*T883;
T909 = getPowerDeriv(T117,(-1),1);
T915 = getPowerDeriv(T342,1-params(3),1);
T919 = getPowerDeriv(T361,1-exp(y(61)),1);
T986 = exp((-((y(33)-T43*T48)/params(13)*(y(33)-T43*T48)/params(13)))/2)/2.506628274631;
T987 = 1/params(13);
T1005 = (-(T987*exp((-(y(33)/params(13)*y(33)/params(13)))/2)/2.506628274631));
T1006 = (1-exp(y(1)))*T1005;
T1018 = exp((-((y(34)-T48*T116)/params(13)*(y(34)-T48*T116)/params(13)))/2)/2.506628274631;
T1036 = (-(T987*exp((-(y(34)/params(13)*y(34)/params(13)))/2)/2.506628274631));
T1037 = (1-exp(y(2)))*T1036;
T1049 = exp((-((y(35)-T43*T48)/params(13)*(y(35)-T43*T48)/params(13)))/2)/2.506628274631;
T1066 = T987*exp((-(y(35)/params(13)*y(35)/params(13)))/2)/2.506628274631;
T1068 = exp(y(1))*(-T1066);
T1080 = exp((-((y(36)-T48*T116)/params(13)*(y(36)-T48*T116)/params(13)))/2)/2.506628274631;
T1097 = T987*exp((-(y(36)/params(13)*y(36)/params(13)))/2)/2.506628274631;
T1099 = exp(y(2))*(-T1097);
T1257 = exp(y(57))*getPowerDeriv(exp(y(57)),1-exp(y(61)),1);
T1259 = T118*T108*T1257;
T1260 = T123*T1259;
T1275 = exp(y(58))*getPowerDeriv(exp(y(58)),1-T32,1);
T1277 = T45*T35*T1275;
T1278 = T52*T1277;
T1339 = (exp(y(61))*exp(y(28))*T116-exp(y(61))*exp(y(61))*exp(y(28)))/(T116*T116);
T1348 = T118*(T112*params(17)*(-exp(y(61)))/(exp(y(61))*exp(y(61)))+T108*T112*(-exp(y(61)))*log(exp(y(57))))+T108*T112*T118*((-exp(y(61)))*log(T117)+(1-exp(y(61)))*T1339/T117);
T1353 = T123*T48*exp(y(61))*2*T116/2;
T1356 = T123*T1348+T119*T1353;
T1448 = (T32*exp(y(27))*T43-T32*T32*exp(y(27)))/(T43*T43);
T1457 = T45*(T39*params(17)*(-T32)/(T32*T32)+T35*T39*(-T32)*log(exp(y(58))))+T35*T39*T45*((-T32)*log(T44)+(1-T32)*T1448/T44);
T1462 = T52*T48*T32*2*T43/2;
T1465 = T52*T1457+T46*T1462;
  g1(1,13)=T10*exp(y(13)-y(27));
  g1(1,15)=exp(y(15));
  g1(1,27)=T10*(-exp(y(13)-y(27)));
  g1(2,14)=T10*exp(y(14)-y(28));
  g1(2,16)=exp(y(16));
  g1(2,28)=T10*(-exp(y(14)-y(28)));
  g1(3,13)=(-(T99*exp(y(25))*params(2)*T90));
  g1(3,85)=(-(T99*exp(y(25))*T90*(-params(2))));
  g1(3,14)=(-T71);
  g1(3,19)=(-(T53*T57*T32*T70));
  g1(3,25)=(-(T91*T99));
  g1(3,27)=(-(T70*T57*T826-exp(y(27)+y(53))*params(15)+T99*exp(y(25))*(1-params(18))*(1-params(2))*T90));
  g1(3,88)=(-(T99*exp(y(25))*T90*(-((1-params(18))*(1-params(2))))));
  g1(3,29)=exp(y(29));
  g1(3,90)=(-(T91*(exp(y(90))+exp(y(53))*(-exp(y(90))))));
  g1(3,92)=(-(T91*exp(y(53))*exp(y(92))));
  g1(3,33)=(-(T70*T53*(-(T986*T987))));
  g1(3,39)=(-(T53*T57*(T32-params(4))*T70));
  g1(3,53)=(-((-(exp(y(27)+y(53))*params(15)))+T91*exp(y(53))*(exp(y(92))-exp(y(90)))));
  g1(3,58)=(-(T70*T57*T1278));
  g1(3,59)=(-(T53*T57*T43*T70));
  g1(3,62)=(-(T70*(T57*T1465+T53*(-(T986*(-(T32*T48/params(13))))))+T53*T57*T70*(T32*y(19)+T32*y(59)+T32*y(39))));
  g1(4,13)=(-T141);
  g1(4,14)=(-(T165*exp(y(26))*params(2)*T156));
  g1(4,86)=(-(T165*exp(y(26))*T156*(-params(2))));
  g1(4,19)=(-(T124*T128*(-exp(y(61)))*T140));
  g1(4,26)=(-T166);
  g1(4,28)=(-(T140*T128*T884-params(15)*exp(y(28)+y(54))+T165*exp(y(26))*(1-params(18))*(1-params(2))*T156));
  g1(4,89)=(-(T165*exp(y(26))*T156*(-((1-params(18))*(1-params(2))))));
  g1(4,30)=exp(y(30));
  g1(4,91)=(-(exp(y(26))*T156*(exp(y(91))+exp(y(54))*(-exp(y(91))))));
  g1(4,93)=(-(exp(y(26))*T156*exp(y(54))*exp(y(93))));
  g1(4,34)=(-(T140*T124*(-(T987*T1018))));
  g1(4,40)=(-(T124*T128*(exp(y(61))-params(4))*T140));
  g1(4,54)=(-((-(params(15)*exp(y(28)+y(54))))+exp(y(26))*T156*exp(y(54))*(exp(y(93))-exp(y(91)))));
  g1(4,57)=(-(T140*T128*T1260));
  g1(4,60)=(-(T124*T128*T116*T140));
  g1(4,61)=(-(T140*(T128*T1356+T124*(-(T1018*(-(T48*exp(y(61))/params(13))))))+T124*T128*T140*(y(19)*(-exp(y(61)))+exp(y(61))*y(60)+exp(y(61))*y(40))));
  g1(5,13)=(-(T184*exp(y(25))*params(2)*T90));
  g1(5,85)=(-(T184*exp(y(25))*T90*(-params(2))));
  g1(5,14)=(-(T70*T53*T173));
  g1(5,19)=(-(T53*T173*T32*T70));
  g1(5,25)=(-(T91*T184));
  g1(5,27)=(-(T70*T173*T826-exp(y(27))*(1-exp(y(55)))*params(16)+T184*exp(y(25))*(1-params(18))*(1-params(2))*T90));
  g1(5,88)=(-(T184*exp(y(25))*T90*(-((1-params(18))*(1-params(2))))));
  g1(5,90)=(-(T91*(exp(y(90))+(1-exp(y(55)))*(-exp(y(90))))));
  g1(5,31)=exp(y(31));
  g1(5,92)=(-(T91*exp(y(92))*(1-exp(y(55)))));
  g1(5,35)=(-(T70*T53*(-(T987*T1049))));
  g1(5,39)=(-(T53*T173*(T32-params(4))*T70));
  g1(5,55)=(-((-(params(16)*exp(y(27))*(-exp(y(55)))))+T91*(exp(y(92))-exp(y(90)))*(-exp(y(55)))));
  g1(5,58)=(-(T70*T173*T1278));
  g1(5,59)=(-(T53*T173*T43*T70));
  g1(5,62)=(-(T53*T173*T70*(T32*y(19)+T32*y(59)+T32*y(39))+T70*(T173*T1465+T53*(-(T1049*(-(T32*T48/params(13))))))));
  g1(6,13)=(-(T140*T124*T192));
  g1(6,14)=(-(T202*exp(y(26))*params(2)*T156));
  g1(6,86)=(-(T202*exp(y(26))*T156*(-params(2))));
  g1(6,19)=(-(T124*T192*(-exp(y(61)))*T140));
  g1(6,26)=(-(exp(y(26))*T156*T202));
  g1(6,28)=(-(T140*T192*T884-params(16)*exp(y(28))*(1-exp(y(56)))+T202*exp(y(26))*(1-params(18))*(1-params(2))*T156));
  g1(6,89)=(-(T202*exp(y(26))*T156*(-((1-params(18))*(1-params(2))))));
  g1(6,91)=(-(exp(y(26))*T156*(exp(y(91))+(1-exp(y(56)))*(-exp(y(91))))));
  g1(6,32)=exp(y(32));
  g1(6,93)=(-(exp(y(26))*T156*exp(y(93))*(1-exp(y(56)))));
  g1(6,36)=(-(T140*T124*(-(T987*T1080))));
  g1(6,40)=(-(T124*T192*(exp(y(61))-params(4))*T140));
  g1(6,56)=(-((-(params(16)*exp(y(28))*(-exp(y(56)))))+exp(y(26))*T156*(exp(y(93))-exp(y(91)))*(-exp(y(56)))));
  g1(6,57)=(-(T140*T192*T1260));
  g1(6,60)=(-(T124*T192*T116*T140));
  g1(6,61)=(-(T124*T192*T140*(y(19)*(-exp(y(61)))+exp(y(61))*y(60)+exp(y(61))*y(40))+T140*(T192*T1356+T124*(-(T1080*(-(T48*exp(y(61))/params(13))))))));
  g1(7,13)=(-((exp(y(92))-exp(y(90)))*exp(y(25))*params(2)*T90));
  g1(7,85)=(-((exp(y(92))-exp(y(90)))*exp(y(25))*T90*(-params(2))));
  g1(7,14)=(-(T46*T212));
  g1(7,19)=(-(T46*T32*T212));
  g1(7,25)=(-(T91*(exp(y(92))-exp(y(90)))));
  g1(7,27)=exp(y(27))*params(15)-(T212*T825+(exp(y(92))-exp(y(90)))*exp(y(25))*(1-params(18))*(1-params(2))*T90);
  g1(7,88)=(-((exp(y(92))-exp(y(90)))*exp(y(25))*T90*(-((1-params(18))*(1-params(2))))));
  g1(7,90)=(-(T91*(-exp(y(90)))));
  g1(7,92)=(-(T91*exp(y(92))));
  g1(7,33)=(-(T46*T43*T212));
  g1(7,39)=(-(T46*(T32-params(4))*T212));
  g1(7,58)=(-(T212*T1277));
  g1(7,59)=(-(T46*T43*T212));
  g1(7,62)=(-(T212*T1457+T46*T212*(T32*y(39)+T32*y(19)+T32*(y(33)+y(59)))));
  g1(8,13)=(-(T119*T223));
  g1(8,14)=(-((exp(y(93))-exp(y(91)))*exp(y(26))*params(2)*T156));
  g1(8,86)=(-((exp(y(93))-exp(y(91)))*exp(y(26))*T156*(-params(2))));
  g1(8,19)=(-(T119*(-exp(y(61)))*T223));
  g1(8,26)=(-(exp(y(26))*T156*(exp(y(93))-exp(y(91)))));
  g1(8,28)=params(15)*exp(y(28))-(T223*T883+(exp(y(93))-exp(y(91)))*exp(y(26))*(1-params(18))*(1-params(2))*T156);
  g1(8,89)=(-((exp(y(93))-exp(y(91)))*exp(y(26))*T156*(-((1-params(18))*(1-params(2))))));
  g1(8,91)=(-(exp(y(26))*T156*(-exp(y(91)))));
  g1(8,93)=(-(exp(y(26))*T156*exp(y(93))));
  g1(8,34)=(-(T119*T116*T223));
  g1(8,40)=(-(T119*(exp(y(61))-params(4))*T223));
  g1(8,57)=(-(T223*T1259));
  g1(8,60)=(-(T119*T116*T223));
  g1(8,61)=(-(T223*T1348+T119*T223*(exp(y(61))*y(40)+y(19)*(-exp(y(61)))+exp(y(61))*(y(34)+y(60)))));
  g1(9,13)=(-((exp(y(92))-exp(y(90)))*exp(y(25))*params(2)*T90));
  g1(9,85)=(-((exp(y(92))-exp(y(90)))*exp(y(25))*T90*(-params(2))));
  g1(9,14)=(-(T46*T234));
  g1(9,19)=(-(T46*T32*T234));
  g1(9,25)=(-(T91*(exp(y(92))-exp(y(90)))));
  g1(9,27)=exp(y(27))*params(16)-((exp(y(92))-exp(y(90)))*exp(y(25))*(1-params(18))*(1-params(2))*T90+T234*T825);
  g1(9,88)=(-((exp(y(92))-exp(y(90)))*exp(y(25))*T90*(-((1-params(18))*(1-params(2))))));
  g1(9,90)=(-(T91*(-exp(y(90)))));
  g1(9,92)=(-(T91*exp(y(92))));
  g1(9,35)=(-(T46*T43*T234));
  g1(9,39)=(-(T46*(T32-params(4))*T234));
  g1(9,58)=(-(T234*T1277));
  g1(9,59)=(-(T46*T43*T234));
  g1(9,62)=(-(T234*T1457+T46*T234*(T32*y(39)+T32*y(19)+T32*(y(59)+y(35)))));
  g1(10,13)=(-(T119*T244));
  g1(10,14)=(-((exp(y(93))-exp(y(91)))*exp(y(26))*params(2)*T156));
  g1(10,86)=(-((exp(y(93))-exp(y(91)))*exp(y(26))*T156*(-params(2))));
  g1(10,19)=(-(T119*(-exp(y(61)))*T244));
  g1(10,26)=(-(exp(y(26))*T156*(exp(y(93))-exp(y(91)))));
  g1(10,28)=exp(y(28))*params(16)-((exp(y(93))-exp(y(91)))*exp(y(26))*(1-params(18))*(1-params(2))*T156+T244*T883);
  g1(10,89)=(-((exp(y(93))-exp(y(91)))*exp(y(26))*T156*(-((1-params(18))*(1-params(2))))));
  g1(10,91)=(-(exp(y(26))*T156*(-exp(y(91)))));
  g1(10,93)=(-(exp(y(26))*T156*exp(y(93))));
  g1(10,36)=(-(T119*T116*T244));
  g1(10,40)=(-(T119*(exp(y(61))-params(4))*T244));
  g1(10,57)=(-(T244*T1259));
  g1(10,60)=(-(T119*T116*T244));
  g1(10,61)=(-(T244*T1348+T119*T244*(exp(y(61))*y(40)+y(19)*(-exp(y(61)))+exp(y(61))*(y(60)+y(36)))));
  g1(11,1)=(-(T258+(1-normcdf(y(33),0,params(13)))*(-exp(y(1)))));
  g1(11,23)=exp(y(23));
  g1(11,33)=(-T1006);
  g1(11,35)=(-T1068);
  g1(12,2)=(-(T271+(1-normcdf(y(34),0,params(13)))*(-exp(y(2)))));
  g1(12,24)=exp(y(24));
  g1(12,34)=(-T1037);
  g1(12,36)=(-T1099);
  g1(13,15)=exp(y(15));
  g1(13,17)=(-exp(y(17)));
  g1(13,1)=(-(params(16)*T258+params(15)*(1-normcdf(y(33),0,params(13)))*(-exp(y(1)))));
  g1(13,33)=(-(params(15)*T1006));
  g1(13,35)=(-(params(16)*T1068));
  g1(14,16)=exp(y(16));
  g1(14,18)=(-exp(y(18)));
  g1(14,2)=(-(params(16)*T271+params(15)*(1-normcdf(y(34),0,params(13)))*(-exp(y(2)))));
  g1(14,34)=(-(params(15)*T1037));
  g1(14,36)=(-(params(16)*T1099));
  g1(15,13)=(-T298);
  g1(15,14)=(-T305);
  g1(15,17)=exp(y(17));
  g1(15,19)=(-T305);
  g1(15,27)=(-(exp(y(13)+(1-params(4))*y(37))*T291*getPowerDeriv(T291,(-1),1)+T304*params(17)*T44*T851));
  g1(15,37)=(-(T292*(1-params(4))*exp(y(13)+(1-params(4))*y(37))));
  g1(15,39)=(-(T300*(1-params(4))*T304));
  g1(15,62)=(-(T304*params(17)*T851*T1448));
  g1(16,13)=(-T323);
  g1(16,14)=(-T315);
  g1(16,18)=exp(y(18));
  g1(16,19)=(-(T317*(-exp(y(13)+(-y(19))+y(40)*(1-params(4))))));
  g1(16,28)=(-(exp(y(14)+(1-params(4))*y(38))*T309*getPowerDeriv(T309,(-1),1)+exp(y(13)+(-y(19))+y(40)*(1-params(4)))*params(17)*T117*T909));
  g1(16,38)=(-(T310*(1-params(4))*exp(y(14)+(1-params(4))*y(38))));
  g1(16,40)=(-(T317*(1-params(4))*exp(y(13)+(-y(19))+y(40)*(1-params(4)))));
  g1(16,61)=(-(exp(y(13)+(-y(19))+y(40)*(1-params(4)))*params(17)*T909*T1339));
  g1(17,27)=(-(params(14)*T332*T857));
  g1(17,37)=(1-params(3))*exp(y(37)*(1-params(3)));
  g1(17,59)=(-(params(14)*T857*params(3)*(-exp(y(27)-y(59)))/(params(3)-1)));
  g1(18,28)=(-(params(14)*T342*T915));
  g1(18,38)=(1-params(3))*exp(y(38)*(1-params(3)));
  g1(18,60)=(-(params(14)*T915*params(3)*(-exp(y(28)-y(60)))/(params(3)-1)));
  g1(19,19)=(1-T32)*exp((1-T32)*(y(19)+y(39)));
  g1(19,27)=(-(exp(y(41))*T39*T350*T861));
  g1(19,39)=(1-T32)*exp((1-T32)*(y(19)+y(39)));
  g1(19,41)=(-T355);
  g1(19,58)=(-(exp(y(41))*T351*T1275));
  g1(19,59)=(-(exp(y(41))*T39*T861*T32*(-exp(y(27)-y(59)))/T43));
  g1(19,62)=exp((1-T32)*(y(19)+y(39)))*(y(19)+y(39))*(-T32)-exp(y(41))*(T351*T39*(-T32)*log(exp(y(58)))+T39*T351*((-T32)*log(T350)+(1-T32)*(T43*T32*exp(y(27)-y(59))-T32*T32*exp(y(27)-y(59)))/(T43*T43)/T350));
  g1(20,19)=exp((1-exp(y(61)))*(y(40)-y(19)))*(-(1-exp(y(61))));
  g1(20,28)=(-(exp(y(42))*T112*T361*T919));
  g1(20,40)=(1-exp(y(61)))*exp((1-exp(y(61)))*(y(40)-y(19)));
  g1(20,42)=(-T366);
  g1(20,57)=(-(exp(y(42))*T362*T1257));
  g1(20,60)=(-(exp(y(42))*T112*T919*exp(y(61))*(-exp(y(28)-y(60)))/T116));
  g1(20,61)=exp((1-exp(y(61)))*(y(40)-y(19)))*(-exp(y(61)))*(y(40)-y(19))-exp(y(42))*(T362*T112*(-exp(y(61)))*log(exp(y(57)))+T112*T362*((-exp(y(61)))*log(T361)+(1-exp(y(61)))*(T116*exp(y(61))*exp(y(28)-y(60))-exp(y(61))*exp(y(61))*exp(y(28)-y(60)))/(T116*T116)/T361));
  g1(21,1)=(-(T173*T52*exp(y(1))+T57*T52*(-exp(y(1)))));
  g1(21,33)=(-(T52*(1-exp(y(1)))*(-(T986*T987))));
  g1(21,35)=(-(T52*exp(y(1))*(-(T987*T1049))));
  g1(21,41)=exp(y(41));
  g1(21,62)=(-(T52*exp(y(1))*(-(T1049*(-(T32*T48/params(13)))))+T173*exp(y(1))*T1462+T52*(1-exp(y(1)))*(-(T986*(-(T32*T48/params(13)))))+T57*(1-exp(y(1)))*T1462));
  g1(22,2)=(-(T192*T123*exp(y(2))+T128*T123*(-exp(y(2)))));
  g1(22,34)=(-(T123*(1-exp(y(2)))*(-(T987*T1018))));
  g1(22,36)=(-(T123*exp(y(2))*(-(T987*T1080))));
  g1(22,42)=exp(y(42));
  g1(22,61)=(-(T123*exp(y(2))*(-(T1080*(-(T48*exp(y(61))/params(13)))))+T192*exp(y(2))*T1353+T123*(1-exp(y(2)))*(-(T1018*(-(T48*exp(y(61))/params(13)))))+T128*(1-exp(y(2)))*T1353));
  g1(23,37)=(-((1-params(4))*exp((1-params(4))*y(37))));
  g1(23,40)=(-(params(17)*(1-params(4))*exp(y(40)*(1-params(4)))));
  g1(24,38)=(-((1-params(4))*exp((1-params(4))*y(38))));
  g1(24,39)=(-(params(17)*(1-params(4))*exp(y(39)*(1-params(4)))));
  g1(25,13)=params(6);
  g1(25,3)=(-params(5));
  g1(25,25)=1;
  g1(25,98)=params(35)/2;
  g1(26,14)=params(6);
  g1(26,4)=(-params(5));
  g1(26,26)=1;
  g1(26,98)=(-(params(35)/2));
  g1(27,13)=(-(exp(y(25))*params(2)*T90));
  g1(27,85)=(-(exp(y(25))*T90*(-params(2))));
  g1(27,25)=(-T91);
  g1(27,27)=(-(exp(y(25))*(1-params(18))*(1-params(2))*T90));
  g1(27,88)=(-(exp(y(25))*T90*(-((1-params(18))*(1-params(2))))));
  g1(27,43)=exp(y(43))*T427+exp(y(43))*params(7)*y(44)*exp(y(43)-y(45));
  g1(27,44)=exp(y(43))*params(7)*exp(y(43)-y(45));
  g1(27,45)=exp(y(43))*params(7)*y(44)*(-exp(y(43)-y(45)));
  g1(28,13)=exp(y(13))-exp(y(13)+(1-params(4))*y(37))*1/params(3);
  g1(28,14)=(-(T35*T304));
  g1(28,17)=(-exp(y(27)+y(17)));
  g1(28,19)=(-(T35*T304));
  g1(28,27)=(-exp(y(27)+y(17)));
  g1(28,37)=(-(1/params(3)*(1-params(4))*exp(y(13)+(1-params(4))*y(37))));
  g1(28,39)=(-(T35*(1-params(4))*T304));
  g1(28,43)=exp(y(43))*y(44);
  g1(28,5)=(-1);
  g1(28,44)=exp(y(43));
  g1(28,62)=(-(T304*params(17)*(-T32)/(T32*T32)));
  g1(29,13)=exp(y(25))*params(2)*T90/T427;
  g1(29,85)=exp(y(25))*T90*(-params(2))/T427;
  g1(29,14)=(-(exp(y(26))*params(2)*T448/T455));
  g1(29,86)=(-(exp(y(26))*T448*(-params(2))/T455));
  g1(29,19)=(-((T449*T455-T449*(-(params(7)*y(44)*(-exp(y(43)-y(46)-y(19))))))/(T455*T455)));
  g1(29,87)=(-(exp(y(26))*(-T448)/T455));
  g1(29,25)=T91/T427;
  g1(29,26)=(-(T449/T455));
  g1(29,27)=exp(y(25))*(1-params(18))*(1-params(2))*T90/T427;
  g1(29,88)=exp(y(25))*T90*(-((1-params(18))*(1-params(2))))/T427;
  g1(29,28)=(-(exp(y(26))*(1-params(18))*(1-params(2))*T448/T455));
  g1(29,89)=(-(exp(y(26))*T448*(-((1-params(18))*(1-params(2))))/T455));
  g1(29,43)=(-(T91*params(7)*y(44)*exp(y(43)-y(45))))/(T427*T427)-(-(T449*(-(params(7)*y(44)*exp(y(43)-y(46)-y(19))))))/(T455*T455);
  g1(29,44)=(-(T91*params(7)*exp(y(43)-y(45))))/(T427*T427)-(-(T449*(-(params(7)*exp(y(43)-y(46)-y(19))))))/(T455*T455);
  g1(29,45)=(-(T91*params(7)*y(44)*(-exp(y(43)-y(45)))))/(T427*T427);
  g1(29,46)=(-((-(T449*(-(params(7)*y(44)*(-exp(y(43)-y(46)-y(19)))))))/(T455*T455)));
  g1(30,22)=1;
  g1(30,45)=(-((-((exp(y(47))-exp(y(48)))*exp(y(45))))/(exp(y(45))*exp(y(45)))));
  g1(30,47)=(-(exp(y(47))/exp(y(45))));
  g1(30,48)=(-((-exp(y(48)))/exp(y(45))));
  g1(31,13)=(-exp(y(13)+(1-params(4))*y(37)));
  g1(31,14)=(-(params(17)*T304));
  g1(31,19)=(-(params(17)*T304));
  g1(31,37)=(-((1-params(4))*exp(y(13)+(1-params(4))*y(37))));
  g1(31,39)=(-(params(17)*(1-params(4))*T304));
  g1(31,45)=exp(y(45));
  g1(32,13)=(-(params(17)*exp(y(13)+(-y(19))+y(40)*(1-params(4)))));
  g1(32,14)=(-exp(y(14)+(1-params(4))*y(38)));
  g1(32,19)=(-(params(17)*(-exp(y(13)+(-y(19))+y(40)*(1-params(4))))));
  g1(32,38)=(-((1-params(4))*exp(y(14)+(1-params(4))*y(38))));
  g1(32,40)=(-(params(17)*(1-params(4))*exp(y(13)+(-y(19))+y(40)*(1-params(4)))));
  g1(32,46)=exp(y(46));
  g1(33,11)=exp(y(11));
  g1(33,37)=exp(y(45)-y(37));
  g1(33,45)=(-exp(y(45)-y(37)));
  g1(34,12)=exp(y(12));
  g1(34,38)=exp(y(46)-y(38));
  g1(34,46)=(-exp(y(46)-y(38)));
  g1(35,14)=(-1);
  g1(35,19)=(-1);
  g1(35,39)=(-(1-params(4)));
  g1(35,47)=1;
  g1(36,13)=(-1);
  g1(36,40)=(-(1-params(4)));
  g1(36,48)=1;
  g1(37,19)=(-T498);
  g1(37,23)=(-(1/T43*T498));
  g1(37,39)=(-T498);
  g1(37,49)=exp(y(49));
  g1(37,58)=T498;
  g1(37,62)=(-(T498*y(23)*(-T32)/(T43*T43)));
  g1(37,99)=(-(params(39)/2));
  g1(38,24)=(-(1/T116*T511));
  g1(38,40)=(-T511);
  g1(38,50)=exp(y(50));
  g1(38,57)=T511;
  g1(38,61)=(-(T511*y(24)*(-exp(y(61)))/(T116*T116)));
  g1(38,99)=params(39)/2;
  g1(39,47)=(-exp(y(47)-y(49)));
  g1(39,49)=exp(y(47)-y(49));
  g1(39,51)=exp(y(51));
  g1(40,48)=(-exp(y(48)-y(50)));
  g1(40,50)=exp(y(48)-y(50));
  g1(40,52)=exp(y(52));
  g1(41,33)=T1005;
  g1(41,53)=(-exp(y(53)));
  g1(42,34)=T1036;
  g1(42,54)=(-exp(y(54)));
  g1(43,35)=T1066;
  g1(43,55)=(-exp(y(55)));
  g1(44,36)=T1097;
  g1(44,56)=(-exp(y(56)));
  g1(45,57)=1;
  g1(45,69)=(-1);
  g1(45,70)=(-0.5);
  g1(45,80)=(-1);
  g1(46,58)=1;
  g1(46,69)=(-1);
  g1(46,70)=0.5;
  g1(47,8)=(-params(8));
  g1(47,75)=1;
  g1(47,94)=(-params(9));
  g1(48,9)=(-params(23));
  g1(48,76)=1;
  g1(48,95)=(-params(22));
  g1(49,59)=1;
  g1(49,75)=(-1);
  g1(49,76)=(-0.5);
  g1(50,60)=1;
  g1(50,75)=(-1);
  g1(50,76)=0.5;
  g1(51,19)=(-(params(3)*params(12)*exp(y(19)*params(12))));
  g1(51,61)=exp(y(61));
  g1(52,19)=(-(params(3)*(-params(12))*exp(y(19)*(-params(12)))));
  g1(52,62)=T32;
  g1(53,20)=1;
  g1(53,49)=1;
  g1(53,50)=(-1);
  g1(54,21)=1;
  g1(54,51)=(-1);
  g1(54,52)=1;
  g1(55,49)=1;
  g1(55,50)=(-1);
  g1(55,63)=1;
  g1(56,45)=(-((-(exp(y(45))*(exp(y(47))+exp(y(48)))))/(exp(y(45))*exp(y(45)))));
  g1(56,47)=(-(exp(y(47))/exp(y(45))));
  g1(56,48)=(-(exp(y(48))/exp(y(45))));
  g1(56,64)=1;
  g1(57,11)=(-((-(exp(y(11))*(exp(y(51))-exp(y(52)))))/(exp(y(11))*exp(y(11)))));
  g1(57,51)=(-(exp(y(51))/exp(y(11))));
  g1(57,52)=(-((-exp(y(52)))/exp(y(11))));
  g1(57,65)=1;
  g1(58,13)=(-1);
  g1(58,14)=1;
  g1(58,21)=(-1);
  g1(58,66)=1;
  g1(58,72)=params(4);
  g1(59,11)=(-1);
  g1(59,12)=1;
  g1(59,21)=(-1);
  g1(59,67)=1;
  g1(59,72)=params(4);
  g1(60,11)=1;
  g1(60,52)=(-1);
  g1(60,68)=1;
  g1(61,6)=(-params(10));
  g1(61,69)=1;
  g1(61,96)=(-(params(11)/params(4)));
  g1(62,7)=(-params(20));
  g1(62,70)=1;
  g1(62,97)=(-(params(21)/params(4)));
  g1(63,11)=1;
  g1(63,12)=(-1);
  g1(63,71)=1;
  g1(64,19)=(-1);
  g1(64,20)=(-1);
  g1(64,72)=1;
  g1(65,11)=(-1);
  g1(65,73)=1;
  g1(66,13)=1;
  g1(66,14)=(-1);
  g1(66,74)=1;
  g1(67,44)=(-(1/exp(y(45))));
  g1(67,45)=(-((-(y(44)*exp(y(45))))/(exp(y(45))*exp(y(45)))));
  g1(67,77)=1;
  g1(68,47)=(-1);
  g1(68,48)=1;
  g1(68,78)=1;
  g1(69,72)=(-1);
  g1(69,79)=1;
  g1(70,74)=(-1);
  g1(70,82)=1;
  g1(71,10)=(-params(38));
  g1(71,80)=1;
  g1(72,11)=(-((-(exp(y(11))*(exp(y(51))+exp(y(52)))))/(exp(y(11))*exp(y(11)))));
  g1(72,51)=(-(exp(y(51))/exp(y(11))));
  g1(72,52)=(-(exp(y(52))/exp(y(11))));
  g1(72,81)=exp(y(81));
  g1(73,71)=(-1);
  g1(73,83)=1;
  g1(74,63)=(-1);
  g1(74,84)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],74,9801);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],74,970299);
end
end
end
end
