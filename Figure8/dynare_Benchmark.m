%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'dynare_Benchmark';
M_.dynare_version = '4.5.6';
oo_.dynare_version = '4.5.6';
options_.dynare_version = '4.5.6';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('dynare_Benchmark.log');
M_.exo_names = 'ezc';
M_.exo_names_tex = 'ezc';
M_.exo_names_long = 'ezc';
M_.exo_names = char(M_.exo_names, 'ezd');
M_.exo_names_tex = char(M_.exo_names_tex, 'ezd');
M_.exo_names_long = char(M_.exo_names_long, 'ezd');
M_.exo_names = char(M_.exo_names, 'exic');
M_.exo_names_tex = char(M_.exo_names_tex, 'exic');
M_.exo_names_long = char(M_.exo_names_long, 'exic');
M_.exo_names = char(M_.exo_names, 'exid');
M_.exo_names_tex = char(M_.exo_names_tex, 'exid');
M_.exo_names_long = char(M_.exo_names_long, 'exid');
M_.exo_names = char(M_.exo_names, 'ebf');
M_.exo_names_tex = char(M_.exo_names_tex, 'ebf');
M_.exo_names_long = char(M_.exo_names_long, 'ebf');
M_.exo_names = char(M_.exo_names, 'etot');
M_.exo_names_tex = char(M_.exo_names_tex, 'etot');
M_.exo_names_long = char(M_.exo_names_long, 'etot');
M_.endo_names = 'YRh';
M_.endo_names_tex = 'YRh';
M_.endo_names_long = 'YRh';
M_.endo_names = char(M_.endo_names, 'YRf');
M_.endo_names_tex = char(M_.endo_names_tex, 'YRf');
M_.endo_names_long = char(M_.endo_names_long, 'YRf');
M_.endo_names = char(M_.endo_names, 'Ch');
M_.endo_names_tex = char(M_.endo_names_tex, 'Ch');
M_.endo_names_long = char(M_.endo_names_long, 'Ch');
M_.endo_names = char(M_.endo_names, 'Cf');
M_.endo_names_tex = char(M_.endo_names_tex, 'Cf');
M_.endo_names_long = char(M_.endo_names_long, 'Cf');
M_.endo_names = char(M_.endo_names, 'Lh');
M_.endo_names_tex = char(M_.endo_names_tex, 'Lh');
M_.endo_names_long = char(M_.endo_names_long, 'Lh');
M_.endo_names = char(M_.endo_names, 'Lf');
M_.endo_names_tex = char(M_.endo_names_tex, 'Lf');
M_.endo_names_long = char(M_.endo_names_long, 'Lf');
M_.endo_names = char(M_.endo_names, 'Lph');
M_.endo_names_tex = char(M_.endo_names_tex, 'Lph');
M_.endo_names_long = char(M_.endo_names_long, 'Lph');
M_.endo_names = char(M_.endo_names, 'Lpf');
M_.endo_names_tex = char(M_.endo_names_tex, 'Lpf');
M_.endo_names_long = char(M_.endo_names_long, 'Lpf');
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'q');
M_.endo_names = char(M_.endo_names, 'TOT');
M_.endo_names_tex = char(M_.endo_names_tex, 'TOT');
M_.endo_names_long = char(M_.endo_names_long, 'TOT');
M_.endo_names = char(M_.endo_names, 'EXIMR');
M_.endo_names_tex = char(M_.endo_names_tex, 'EXIMR');
M_.endo_names_long = char(M_.endo_names_long, 'EXIMR');
M_.endo_names = char(M_.endo_names, 'NXY');
M_.endo_names_tex = char(M_.endo_names_tex, 'NXY');
M_.endo_names_long = char(M_.endo_names_long, 'NXY');
M_.endo_names = char(M_.endo_names, 'Nh');
M_.endo_names_tex = char(M_.endo_names_tex, 'Nh');
M_.endo_names_long = char(M_.endo_names_long, 'Nh');
M_.endo_names = char(M_.endo_names, 'Nf');
M_.endo_names_tex = char(M_.endo_names_tex, 'Nf');
M_.endo_names_long = char(M_.endo_names_long, 'Nf');
M_.endo_names = char(M_.endo_names, 'beth');
M_.endo_names_tex = char(M_.endo_names_tex, 'beth');
M_.endo_names_long = char(M_.endo_names_long, 'beth');
M_.endo_names = char(M_.endo_names, 'betf');
M_.endo_names_tex = char(M_.endo_names_tex, 'betf');
M_.endo_names_long = char(M_.endo_names_long, 'betf');
M_.endo_names = char(M_.endo_names, 'Wh');
M_.endo_names_tex = char(M_.endo_names_tex, 'Wh');
M_.endo_names_long = char(M_.endo_names_long, 'Wh');
M_.endo_names = char(M_.endo_names, 'Wf');
M_.endo_names_tex = char(M_.endo_names_tex, 'Wf');
M_.endo_names_long = char(M_.endo_names_long, 'Wf');
M_.endo_names = char(M_.endo_names, 'EV0h');
M_.endo_names_tex = char(M_.endo_names_tex, 'EV0h');
M_.endo_names_long = char(M_.endo_names_long, 'EV0h');
M_.endo_names = char(M_.endo_names, 'EV0f');
M_.endo_names_tex = char(M_.endo_names_tex, 'EV0f');
M_.endo_names_long = char(M_.endo_names_long, 'EV0f');
M_.endo_names = char(M_.endo_names, 'EV1h');
M_.endo_names_tex = char(M_.endo_names_tex, 'EV1h');
M_.endo_names_long = char(M_.endo_names_long, 'EV1h');
M_.endo_names = char(M_.endo_names, 'EV1f');
M_.endo_names_tex = char(M_.endo_names_tex, 'EV1f');
M_.endo_names_long = char(M_.endo_names_long, 'EV1f');
M_.endo_names = char(M_.endo_names, 'eta0h');
M_.endo_names_tex = char(M_.endo_names_tex, 'eta0h');
M_.endo_names_long = char(M_.endo_names_long, 'eta0h');
M_.endo_names = char(M_.endo_names, 'eta0f');
M_.endo_names_tex = char(M_.endo_names_tex, 'eta0f');
M_.endo_names_long = char(M_.endo_names_long, 'eta0f');
M_.endo_names = char(M_.endo_names, 'eta1h');
M_.endo_names_tex = char(M_.endo_names_tex, 'eta1h');
M_.endo_names_long = char(M_.endo_names_long, 'eta1h');
M_.endo_names = char(M_.endo_names, 'eta1f');
M_.endo_names_tex = char(M_.endo_names_tex, 'eta1f');
M_.endo_names_long = char(M_.endo_names_long, 'eta1f');
M_.endo_names = char(M_.endo_names, 'Ph');
M_.endo_names_tex = char(M_.endo_names_tex, 'Ph');
M_.endo_names_long = char(M_.endo_names_long, 'Ph');
M_.endo_names = char(M_.endo_names, 'Pfs');
M_.endo_names_tex = char(M_.endo_names_tex, 'Pfs');
M_.endo_names_long = char(M_.endo_names_long, 'Pfs');
M_.endo_names = char(M_.endo_names, 'Phs');
M_.endo_names_tex = char(M_.endo_names_tex, 'Phs');
M_.endo_names_long = char(M_.endo_names_long, 'Phs');
M_.endo_names = char(M_.endo_names, 'Pf');
M_.endo_names_tex = char(M_.endo_names_tex, 'Pf');
M_.endo_names_long = char(M_.endo_names_long, 'Pf');
M_.endo_names = char(M_.endo_names, 'PsiXh');
M_.endo_names_tex = char(M_.endo_names_tex, 'PsiXh');
M_.endo_names_long = char(M_.endo_names_long, 'PsiXh');
M_.endo_names = char(M_.endo_names, 'PsiXf');
M_.endo_names_tex = char(M_.endo_names_tex, 'PsiXf');
M_.endo_names_long = char(M_.endo_names_long, 'PsiXf');
M_.endo_names = char(M_.endo_names, 'V');
M_.endo_names_tex = char(M_.endo_names_tex, 'V');
M_.endo_names_long = char(M_.endo_names_long, 'V');
M_.endo_names = char(M_.endo_names, 'B');
M_.endo_names_tex = char(M_.endo_names_tex, 'B');
M_.endo_names_long = char(M_.endo_names_long, 'B');
M_.endo_names = char(M_.endo_names, 'YNh');
M_.endo_names_tex = char(M_.endo_names_tex, 'YNh');
M_.endo_names_long = char(M_.endo_names_long, 'YNh');
M_.endo_names = char(M_.endo_names, 'YNf');
M_.endo_names_tex = char(M_.endo_names_tex, 'YNf');
M_.endo_names_long = char(M_.endo_names_long, 'YNf');
M_.endo_names = char(M_.endo_names, 'EXN');
M_.endo_names_tex = char(M_.endo_names_tex, 'EXN');
M_.endo_names_long = char(M_.endo_names_long, 'EXN');
M_.endo_names = char(M_.endo_names, 'IMN');
M_.endo_names_tex = char(M_.endo_names_tex, 'IMN');
M_.endo_names_long = char(M_.endo_names_long, 'IMN');
M_.endo_names = char(M_.endo_names, 'Px');
M_.endo_names_tex = char(M_.endo_names_tex, 'Px');
M_.endo_names_long = char(M_.endo_names_long, 'Px');
M_.endo_names = char(M_.endo_names, 'Pm');
M_.endo_names_tex = char(M_.endo_names_tex, 'Pm');
M_.endo_names_long = char(M_.endo_names_long, 'Pm');
M_.endo_names = char(M_.endo_names, 'EXR');
M_.endo_names_tex = char(M_.endo_names_tex, 'EXR');
M_.endo_names_long = char(M_.endo_names_long, 'EXR');
M_.endo_names = char(M_.endo_names, 'IMR');
M_.endo_names_tex = char(M_.endo_names_tex, 'IMR');
M_.endo_names_long = char(M_.endo_names_long, 'IMR');
M_.endo_names = char(M_.endo_names, 'n0h');
M_.endo_names_tex = char(M_.endo_names_tex, 'n0h');
M_.endo_names_long = char(M_.endo_names_long, 'n0h');
M_.endo_names = char(M_.endo_names, 'n0f');
M_.endo_names_tex = char(M_.endo_names_tex, 'n0f');
M_.endo_names_long = char(M_.endo_names_long, 'n0f');
M_.endo_names = char(M_.endo_names, 'n1h');
M_.endo_names_tex = char(M_.endo_names_tex, 'n1h');
M_.endo_names_long = char(M_.endo_names_long, 'n1h');
M_.endo_names = char(M_.endo_names, 'n1f');
M_.endo_names_tex = char(M_.endo_names_tex, 'n1f');
M_.endo_names_long = char(M_.endo_names_long, 'n1f');
M_.endo_names = char(M_.endo_names, 'xih');
M_.endo_names_tex = char(M_.endo_names_tex, 'xih');
M_.endo_names_long = char(M_.endo_names_long, 'xih');
M_.endo_names = char(M_.endo_names, 'xif');
M_.endo_names_tex = char(M_.endo_names_tex, 'xif');
M_.endo_names_long = char(M_.endo_names_long, 'xif');
M_.endo_names = char(M_.endo_names, 'Zh');
M_.endo_names_tex = char(M_.endo_names_tex, 'Zh');
M_.endo_names_long = char(M_.endo_names_long, 'Zh');
M_.endo_names = char(M_.endo_names, 'Zf');
M_.endo_names_tex = char(M_.endo_names_tex, 'Zf');
M_.endo_names_long = char(M_.endo_names_long, 'Zf');
M_.endo_names = char(M_.endo_names, 'thh');
M_.endo_names_tex = char(M_.endo_names_tex, 'thh');
M_.endo_names_long = char(M_.endo_names_long, 'thh');
M_.endo_names = char(M_.endo_names, 'thf');
M_.endo_names_tex = char(M_.endo_names_tex, 'thf');
M_.endo_names_long = char(M_.endo_names_long, 'thf');
M_.endo_names = char(M_.endo_names, 'TOTa');
M_.endo_names_tex = char(M_.endo_names_tex, 'TOTa');
M_.endo_names_long = char(M_.endo_names_long, 'TOTa');
M_.endo_names = char(M_.endo_names, 'tshare');
M_.endo_names_tex = char(M_.endo_names_tex, 'tshare');
M_.endo_names_long = char(M_.endo_names_long, 'tshare');
M_.endo_names = char(M_.endo_names, 'NXYR');
M_.endo_names_tex = char(M_.endo_names_tex, 'NXYR');
M_.endo_names_long = char(M_.endo_names_long, 'NXYR');
M_.endo_names = char(M_.endo_names, 'WEDGEDD');
M_.endo_names_tex = char(M_.endo_names_tex, 'WEDGEDD');
M_.endo_names_long = char(M_.endo_names_long, 'WEDGEDD');
M_.endo_names = char(M_.endo_names, 'WEDGEYY');
M_.endo_names_tex = char(M_.endo_names_tex, 'WEDGEYY');
M_.endo_names_long = char(M_.endo_names_long, 'WEDGEYY');
M_.endo_names = char(M_.endo_names, 'MD');
M_.endo_names_tex = char(M_.endo_names_tex, 'MD');
M_.endo_names_long = char(M_.endo_names_long, 'MD');
M_.endo_names = char(M_.endo_names, 'xic');
M_.endo_names_tex = char(M_.endo_names_tex, 'xic');
M_.endo_names_long = char(M_.endo_names_long, 'xic');
M_.endo_names = char(M_.endo_names, 'xid');
M_.endo_names_tex = char(M_.endo_names_tex, 'xid');
M_.endo_names_long = char(M_.endo_names_long, 'xid');
M_.endo_names = char(M_.endo_names, 'YSY');
M_.endo_names_tex = char(M_.endo_names_tex, 'YSY');
M_.endo_names_long = char(M_.endo_names_long, 'YSY');
M_.endo_names = char(M_.endo_names, 'totq');
M_.endo_names_tex = char(M_.endo_names_tex, 'totq');
M_.endo_names_long = char(M_.endo_names_long, 'totq');
M_.endo_names = char(M_.endo_names, 'ipdetrend');
M_.endo_names_tex = char(M_.endo_names_tex, 'ipdetrend');
M_.endo_names_long = char(M_.endo_names_long, 'ipdetrend');
M_.endo_names = char(M_.endo_names, 'DSD_ADV');
M_.endo_names_tex = char(M_.endo_names_tex, 'DSD\_ADV');
M_.endo_names_long = char(M_.endo_names_long, 'DSD_ADV');
M_.endo_names = char(M_.endo_names, 'Zc');
M_.endo_names_tex = char(M_.endo_names_tex, 'Zc');
M_.endo_names_long = char(M_.endo_names_long, 'Zc');
M_.endo_names = char(M_.endo_names, 'Zd');
M_.endo_names_tex = char(M_.endo_names_tex, 'Zd');
M_.endo_names_long = char(M_.endo_names_long, 'Zd');
M_.endo_names = char(M_.endo_names, 'By');
M_.endo_names_tex = char(M_.endo_names_tex, 'By');
M_.endo_names_long = char(M_.endo_names_long, 'By');
M_.endo_names = char(M_.endo_names, 'EXIMN');
M_.endo_names_tex = char(M_.endo_names_tex, 'EXIMN');
M_.endo_names_long = char(M_.endo_names_long, 'EXIMN');
M_.endo_names = char(M_.endo_names, 'mtotq');
M_.endo_names_tex = char(M_.endo_names_tex, 'mtotq');
M_.endo_names_long = char(M_.endo_names_long, 'mtotq');
M_.endo_names = char(M_.endo_names, 'xia');
M_.endo_names_tex = char(M_.endo_names_tex, 'xia');
M_.endo_names_long = char(M_.endo_names_long, 'xia');
M_.endo_names = char(M_.endo_names, 'XMY');
M_.endo_names_tex = char(M_.endo_names_tex, 'XMY');
M_.endo_names_long = char(M_.endo_names_long, 'XMY');
M_.endo_names = char(M_.endo_names, 'mDSD_ADV');
M_.endo_names_tex = char(M_.endo_names_tex, 'mDSD\_ADV');
M_.endo_names_long = char(M_.endo_names_long, 'mDSD_ADV');
M_.endo_names = char(M_.endo_names, 'mysy');
M_.endo_names_tex = char(M_.endo_names_tex, 'mysy');
M_.endo_names_long = char(M_.endo_names_long, 'mysy');
M_.endo_names = char(M_.endo_names, 'mtot');
M_.endo_names_tex = char(M_.endo_names_tex, 'mtot');
M_.endo_names_long = char(M_.endo_names_long, 'mtot');
M_.endo_partitions = struct();
M_.param_names = 'bet';
M_.param_names_tex = 'bet';
M_.param_names_long = 'bet';
M_.param_names = char(M_.param_names, 'sig');
M_.param_names_tex = char(M_.param_names_tex, 'sig');
M_.param_names_long = char(M_.param_names_long, 'sig');
M_.param_names = char(M_.param_names, 'th');
M_.param_names_tex = char(M_.param_names_tex, 'th');
M_.param_names_long = char(M_.param_names_long, 'th');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'rhob');
M_.param_names_tex = char(M_.param_names_tex, 'rhob');
M_.param_names_long = char(M_.param_names_long, 'rhob');
M_.param_names = char(M_.param_names, 'fi');
M_.param_names_tex = char(M_.param_names_tex, 'fi');
M_.param_names_long = char(M_.param_names_long, 'fi');
M_.param_names = char(M_.param_names, 'xib');
M_.param_names_tex = char(M_.param_names_tex, 'xib');
M_.param_names_long = char(M_.param_names_long, 'xib');
M_.param_names = char(M_.param_names, 'rhozc');
M_.param_names_tex = char(M_.param_names_tex, 'rhozc');
M_.param_names_long = char(M_.param_names_long, 'rhozc');
M_.param_names = char(M_.param_names, 'sdzc');
M_.param_names_tex = char(M_.param_names_tex, 'sdzc');
M_.param_names_long = char(M_.param_names_long, 'sdzc');
M_.param_names = char(M_.param_names, 'rhoxic');
M_.param_names_tex = char(M_.param_names_tex, 'rhoxic');
M_.param_names_long = char(M_.param_names_long, 'rhoxic');
M_.param_names = char(M_.param_names, 'sdci');
M_.param_names_tex = char(M_.param_names_tex, 'sdci');
M_.param_names_long = char(M_.param_names_long, 'sdci');
M_.param_names = char(M_.param_names, 'zeta');
M_.param_names_tex = char(M_.param_names_tex, 'zeta');
M_.param_names_long = char(M_.param_names_long, 'zeta');
M_.param_names = char(M_.param_names, 'sdeta');
M_.param_names_tex = char(M_.param_names_tex, 'sdeta');
M_.param_names_long = char(M_.param_names_long, 'sdeta');
M_.param_names = char(M_.param_names, 'PsiT');
M_.param_names_tex = char(M_.param_names_tex, 'PsiT');
M_.param_names_long = char(M_.param_names_long, 'PsiT');
M_.param_names = char(M_.param_names, 'tau0');
M_.param_names_tex = char(M_.param_names_tex, 'tau0');
M_.param_names_long = char(M_.param_names_long, 'tau0');
M_.param_names = char(M_.param_names, 'tau1');
M_.param_names_tex = char(M_.param_names_tex, 'tau1');
M_.param_names_long = char(M_.param_names_long, 'tau1');
M_.param_names = char(M_.param_names, 'a2');
M_.param_names_tex = char(M_.param_names_tex, 'a2');
M_.param_names_long = char(M_.param_names_long, 'a2');
M_.param_names = char(M_.param_names, 'gam');
M_.param_names_tex = char(M_.param_names_tex, 'gam');
M_.param_names_long = char(M_.param_names_long, 'gam');
M_.param_names = char(M_.param_names, 'Css');
M_.param_names_tex = char(M_.param_names_tex, 'Css');
M_.param_names_long = char(M_.param_names_long, 'Css');
M_.param_names = char(M_.param_names, 'rhoxid');
M_.param_names_tex = char(M_.param_names_tex, 'rhoxid');
M_.param_names_long = char(M_.param_names_long, 'rhoxid');
M_.param_names = char(M_.param_names, 'sddi');
M_.param_names_tex = char(M_.param_names_tex, 'sddi');
M_.param_names_long = char(M_.param_names_long, 'sddi');
M_.param_names = char(M_.param_names, 'sdzd');
M_.param_names_tex = char(M_.param_names_tex, 'sdzd');
M_.param_names_long = char(M_.param_names_long, 'sdzd');
M_.param_names = char(M_.param_names, 'rhozd');
M_.param_names_tex = char(M_.param_names_tex, 'rhozd');
M_.param_names_long = char(M_.param_names_long, 'rhozd');
M_.param_names = char(M_.param_names, 'Yss');
M_.param_names_tex = char(M_.param_names_tex, 'Yss');
M_.param_names_long = char(M_.param_names_long, 'Yss');
M_.param_names = char(M_.param_names, 'Psi0');
M_.param_names_tex = char(M_.param_names_tex, 'Psi0');
M_.param_names_long = char(M_.param_names_long, 'Psi0');
M_.param_names = char(M_.param_names, 'Psi1');
M_.param_names_tex = char(M_.param_names_tex, 'Psi1');
M_.param_names_long = char(M_.param_names_long, 'Psi1');
M_.param_names = char(M_.param_names, 'PsiX');
M_.param_names_tex = char(M_.param_names_tex, 'PsiX');
M_.param_names_long = char(M_.param_names_long, 'PsiX');
M_.param_names = char(M_.param_names, 'n1');
M_.param_names_tex = char(M_.param_names_tex, 'n1');
M_.param_names_long = char(M_.param_names_long, 'n1');
M_.param_names = char(M_.param_names, 'n0');
M_.param_names_tex = char(M_.param_names_tex, 'n0');
M_.param_names_long = char(M_.param_names_long, 'n0');
M_.param_names = char(M_.param_names, 'N');
M_.param_names_tex = char(M_.param_names_tex, 'N');
M_.param_names_long = char(M_.param_names_long, 'N');
M_.param_names = char(M_.param_names, 'eta0');
M_.param_names_tex = char(M_.param_names_tex, 'eta0');
M_.param_names_long = char(M_.param_names_long, 'eta0');
M_.param_names = char(M_.param_names, 'eta1');
M_.param_names_tex = char(M_.param_names_tex, 'eta1');
M_.param_names_long = char(M_.param_names_long, 'eta1');
M_.param_names = char(M_.param_names, 'Tr');
M_.param_names_tex = char(M_.param_names_tex, 'Tr');
M_.param_names_long = char(M_.param_names_long, 'Tr');
M_.param_names = char(M_.param_names, 'L');
M_.param_names_tex = char(M_.param_names_tex, 'L');
M_.param_names_long = char(M_.param_names_long, 'L');
M_.param_names = char(M_.param_names, 'sdb');
M_.param_names_tex = char(M_.param_names_tex, 'sdb');
M_.param_names_long = char(M_.param_names_long, 'sdb');
M_.param_names = char(M_.param_names, 'atotq');
M_.param_names_tex = char(M_.param_names_tex, 'atotq');
M_.param_names_long = char(M_.param_names_long, 'atotq');
M_.param_names = char(M_.param_names, 'aysy');
M_.param_names_tex = char(M_.param_names_tex, 'aysy');
M_.param_names_long = char(M_.param_names_long, 'aysy');
M_.param_names = char(M_.param_names, 'rhoxia');
M_.param_names_tex = char(M_.param_names_tex, 'rhoxia');
M_.param_names_long = char(M_.param_names_long, 'rhoxia');
M_.param_names = char(M_.param_names, 'sdtot');
M_.param_names_tex = char(M_.param_names_tex, 'sdtot');
M_.param_names_long = char(M_.param_names_long, 'sdtot');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 6;
M_.endo_nbr = 74;
M_.param_nbr = 39;
M_.orig_endo_nbr = 74;
M_.aux_vars = [];
options_.varobs = cell(1);
options_.varobs(1)  = {'EXIMR'};
options_.varobs(2)  = {'XMY'};
options_.varobs(3)  = {'mysy'};
options_.varobs(4)  = {'ipdetrend'};
options_.varobs(5)  = {'mtotq'};
options_.varobs(6)  = {'mtot'};
options_.varobs_id = [ 11 71 73 63 69 74  ];
M_.Sigma_e = zeros(6, 6);
M_.Correlation_matrix = eye(6, 6);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 0;
erase_compiled_function('dynare_Benchmark_static');
erase_compiled_function('dynare_Benchmark_dynamic');
M_.orig_eq_nbr = 74;
M_.eq_nbr = 74;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 0 11 0;
 0 12 0;
 0 13 85;
 0 14 86;
 0 15 0;
 0 16 0;
 0 17 0;
 0 18 0;
 0 19 87;
 0 20 0;
 0 21 0;
 0 22 0;
 1 23 0;
 2 24 0;
 3 25 0;
 4 26 0;
 0 27 88;
 0 28 89;
 0 29 90;
 0 30 91;
 0 31 92;
 0 32 93;
 0 33 0;
 0 34 0;
 0 35 0;
 0 36 0;
 0 37 0;
 0 38 0;
 0 39 0;
 0 40 0;
 0 41 0;
 0 42 0;
 0 43 0;
 5 44 0;
 0 45 0;
 0 46 0;
 0 47 0;
 0 48 0;
 0 49 0;
 0 50 0;
 0 51 0;
 0 52 0;
 0 53 0;
 0 54 0;
 0 55 0;
 0 56 0;
 0 57 0;
 0 58 0;
 0 59 0;
 0 60 0;
 0 61 0;
 0 62 0;
 0 63 0;
 0 64 0;
 0 65 0;
 0 66 0;
 0 67 0;
 0 68 0;
 6 69 0;
 7 70 0;
 0 71 0;
 0 72 0;
 0 73 0;
 0 74 0;
 8 75 0;
 9 76 0;
 0 77 0;
 0 78 0;
 0 79 0;
 10 80 0;
 0 81 0;
 0 82 0;
 0 83 0;
 0 84 0;]';
M_.nstatic = 55;
M_.nfwrd   = 9;
M_.npred   = 10;
M_.nboth   = 0;
M_.nsfwrd   = 9;
M_.nspred   = 10;
M_.ndynamic   = 19;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:6];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(74, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(6, 1);
M_.params = NaN(39, 1);
M_.NNZDerivatives = [383; 1521; -1];
clc;
close all;
load xss;
load xpar;
M_.params( 1 ) = xpar(1);
bet = M_.params( 1 );
M_.params( 2 ) = xpar(2);
sig = M_.params( 2 );
M_.params( 3 ) = xpar(3);
th = M_.params( 3 );
M_.params( 4 ) = xpar(4);
rho = M_.params( 4 );
M_.params( 5 ) = xpar(5);
rhob = M_.params( 5 );
M_.params( 6 ) = xpar(6);
fi = M_.params( 6 );
M_.params( 7 ) = xpar(7);
xib = M_.params( 7 );
M_.params( 8 ) = xpar(8);
rhozc = M_.params( 8 );
M_.params( 9 ) = xpar(9);
sdzc = M_.params( 9 );
M_.params( 10 ) = xpar(10);
rhoxic = M_.params( 10 );
M_.params( 11 ) = xpar(11);
sdci = M_.params( 11 );
M_.params( 12 ) = xpar(12);
zeta = M_.params( 12 );
M_.params( 13 ) = xpar(13);
sdeta = M_.params( 13 );
M_.params( 14 ) = xpar(14);
PsiT = M_.params( 14 );
M_.params( 15 ) = xpar(15);
tau0 = M_.params( 15 );
M_.params( 16 ) = xpar(16);
tau1 = M_.params( 16 );
M_.params( 17 ) = xpar(17);
a2 = M_.params( 17 );
M_.params( 18 ) = xpar(18);
gam = M_.params( 18 );
M_.params( 19 ) = xpar(19);
Css = M_.params( 19 );
M_.params( 20 ) = xpar(20);
rhoxid = M_.params( 20 );
M_.params( 21 ) = xpar(21);
sddi = M_.params( 21 );
M_.params( 22 ) = xpar(22);
sdzd = M_.params( 22 );
M_.params( 23 ) = xpar(23);
rhozd = M_.params( 23 );
M_.params( 24 ) = xpar(24);
Yss = M_.params( 24 );
M_.params( 25 ) = xpar(25);
Psi0 = M_.params( 25 );
M_.params( 26 ) = xpar(26);
Psi1 = M_.params( 26 );
M_.params( 27 ) = xpar(27);
PsiX = M_.params( 27 );
M_.params( 28 ) = xpar(28);
n1 = M_.params( 28 );
M_.params( 29 ) = xpar(29);
n0 = M_.params( 29 );
M_.params( 30 ) = xpar(30);
N = M_.params( 30 );
M_.params( 31 ) = xpar(31);
eta0 = M_.params( 31 );
M_.params( 32 ) = xpar(32);
eta1 = M_.params( 32 );
M_.params( 33 ) = xpar(33);
Tr = M_.params( 33 );
M_.params( 34 ) = xpar(34);
L = M_.params( 34 );
M_.params( 35 ) = xpar(35);
sdb = M_.params( 35 );
M_.params( 36 ) = xpar(36);
atotq = M_.params( 36 );
M_.params( 37 ) = xpar(37);
aysy = M_.params( 37 );
M_.params( 38 ) = xpar(38);
rhoxia = M_.params( 38 );
M_.params( 39 ) = xpar(39);
sdtot = M_.params( 39 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 17 ) = xss(1);
oo_.steady_state( 18 ) = xss(2);
oo_.steady_state( 19 ) = xss(3);
oo_.steady_state( 20 ) = xss(4);
oo_.steady_state( 21 ) = xss(5);
oo_.steady_state( 22 ) = xss(6);
oo_.steady_state( 23 ) = xss(7);
oo_.steady_state( 24 ) = xss(8);
oo_.steady_state( 25 ) = xss(9);
oo_.steady_state( 26 ) = xss(10);
oo_.steady_state( 13 ) = xss(11);
oo_.steady_state( 14 ) = xss(12);
oo_.steady_state( 5 ) = xss(13);
oo_.steady_state( 6 ) = xss(14);
oo_.steady_state( 7 ) = xss(15);
oo_.steady_state( 8 ) = xss(16);
oo_.steady_state( 27 ) = xss(17);
oo_.steady_state( 28 ) = xss(18);
oo_.steady_state( 29 ) = xss(19);
oo_.steady_state( 30 ) = xss(20);
oo_.steady_state( 31 ) = xss(21);
oo_.steady_state( 32 ) = xss(22);
oo_.steady_state( 3 ) = xss(23);
oo_.steady_state( 4 ) = xss(24);
oo_.steady_state( 15 ) = xss(25);
oo_.steady_state( 16 ) = xss(26);
oo_.steady_state( 33 ) = xss(27);
oo_.steady_state( 34 ) = xss(28);
oo_.steady_state( 9 ) = xss(29);
oo_.steady_state( 12 ) = xss(30);
oo_.steady_state( 35 ) = xss(31);
oo_.steady_state( 36 ) = xss(32);
oo_.steady_state( 1 ) = xss(33);
oo_.steady_state( 2 ) = xss(34);
oo_.steady_state( 37 ) = xss(35);
oo_.steady_state( 38 ) = xss(36);
oo_.steady_state( 39 ) = xss(37);
oo_.steady_state( 40 ) = xss(38);
oo_.steady_state( 41 ) = xss(39);
oo_.steady_state( 42 ) = xss(40);
oo_.steady_state( 43 ) = xss(41);
oo_.steady_state( 44 ) = xss(42);
oo_.steady_state( 45 ) = xss(43);
oo_.steady_state( 46 ) = xss(44);
oo_.steady_state( 47 ) = xss(45);
oo_.steady_state( 48 ) = xss(46);
oo_.steady_state( 49 ) = xss(47);
oo_.steady_state( 50 ) = xss(48);
oo_.steady_state( 51 ) = xss(49);
oo_.steady_state( 52 ) = xss(50);
oo_.steady_state( 59 ) = xss(51);
oo_.steady_state( 60 ) = xss(52);
oo_.steady_state( 10 ) = xss(53);
oo_.steady_state( 11 ) = xss(54);
oo_.steady_state( 53 ) = xss(55);
oo_.steady_state( 54 ) = xss(56);
oo_.steady_state( 55 ) = xss(57);
oo_.steady_state( 56 ) = xss(58);
oo_.steady_state( 57 ) = xss(59);
oo_.steady_state( 58 ) = xss(60);
oo_.steady_state( 61 ) = xss(61);
oo_.steady_state( 62 ) = xss(62);
oo_.steady_state( 63 ) = xss(63);
oo_.steady_state( 64 ) = xss(64);
oo_.steady_state( 65 ) = xss(65);
oo_.steady_state( 66 ) = xss(66);
oo_.steady_state( 67 ) = xss(67);
oo_.steady_state( 68 ) = xss(68);
oo_.steady_state( 69 ) = xss(69);
oo_.steady_state( 70 ) = xss(70);
oo_.steady_state( 71 ) = xss(71);
oo_.steady_state( 72 ) = xss(72);
oo_.exo_steady_state( 5 ) = 0;
oo_.steady_state( 73 ) = xss(73);
oo_.steady_state( 74 ) = xss(74);
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 1;
M_.Sigma_e(3, 3) = 1;
M_.Sigma_e(4, 4) = 1;
M_.Sigma_e(5, 5) = 1;
M_.Sigma_e(6, 6) = 1;
M_.Sigma_e(1, 2) = 0.0;
M_.Sigma_e(2, 1) = M_.Sigma_e(1, 2);
M_.Sigma_e(3, 4) = 0.0;
M_.Sigma_e(4, 3) = M_.Sigma_e(3, 4);
M_.sigma_e_is_diagonal = 0;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
resid(1);
estim_params_.var_exo = [];
estim_params_.var_endo = [];
estim_params_.corrx = [];
estim_params_.corrn = [];
estim_params_.param_vals = [];
estim_params_.param_vals = [estim_params_.param_vals; 8, 0.98, (-Inf), Inf, 5, NaN, NaN, 0.9, 1, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 23, 0.95, (-Inf), Inf, 5, NaN, NaN, 0.9, 1, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 10, 0.98, (-Inf), Inf, 5, NaN, NaN, 0.9, 1, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 20, 0.95, (-Inf), Inf, 5, NaN, NaN, 0.9, 1, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 5, 0.95, (-Inf), Inf, 5, NaN, NaN, 0.9, 1, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 11, 0.01, (-Inf), Inf, 4, 0.01, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 21, 0.01, (-Inf), Inf, 4, 0.01, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 9, 0.01, (-Inf), Inf, 4, 0.01, 0.02, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 22, 0.007, (-Inf), Inf, 4, 0.01, 0.02, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 35, 0.001, (-Inf), Inf, 4, 0.001, 0.02, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 39, 0.01, (-Inf), Inf, 4, 0.01, 0.05, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 12, 0.5, (-Inf), Inf, 3, 0.5, 0.1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 6, 0.01, (-Inf), Inf, 4, 0.01, 0.1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 4, 3.3, (-Inf), Inf, 5, NaN, NaN, 1, 3.99, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 2, 6, (-Inf), Inf, 5, NaN, NaN, 1, 8, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 36, 0, (-Inf), Inf, 3, 0, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 37, 0, (-Inf), Inf, 3, 0, 0.2, NaN, NaN, NaN ];
options_.use_calibration_initialization = 1;
options_.mh_jscale = 0.3;
options_.mh_nblck = 5;
options_.mh_replic = 200000;
options_.mode_compute = 6;
options_.plot_priors = 1;
options_.datafile = 'Data02_world';
options_.nobs = 139;
options_.order = 1;
var_list_ = char();
oo_recursive_=dynare_estimation(var_list_);
options_.parameter_set = 'posterior_mode';
var_list_ = char();
[oo_,M_]= shock_decomposition(M_,oo_,options_,var_list_,bayestopt_,estim_params_);
save('dynare_Benchmark_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('dynare_Benchmark_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('dynare_Benchmark_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('dynare_Benchmark_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('dynare_Benchmark_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('dynare_Benchmark_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('dynare_Benchmark_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
