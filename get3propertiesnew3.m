function [sigma_m,ps_m,vl_m] = get3propertiesnew3(T_guess)
% This function, which is called by the host function PlotNucTempPseudo,
% calculates three key properties for the surrogate fuel mixture. These
% three properties are surface tension (sigma_m), saturation pressure
% (ps_m), and liquid-molar volume (vl_m). These are integral to calculating
% the nucleation temperatures of the surrogate mixture, which then can be
% used to predict its superheat limit. 





%%%%%%%%%%%%%%%%%%%%%%%%% Set Up Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Get mixture critical temperature and pressure from another program in
% this bundle.
[Tcm,Pcm] = MixtureCriticalTP(); %critical temperature of surrogate mixture

% Pre-allocate key constants
sigma0 = 0;
mu = 0;
s = 0;
beta = 0;
gamma = 0;
m = 0;
A = 0;
B = 0;
C = 0;
D = 0;
Zra = 0;

% Old code to be turned on/off for modification and testing purposes
%Tc = 616.2;   Tb = 411.5;   Pc = 35.10000;  Zra = 0.2592;  %p-xylene
%Tc = 748.4;   Tb = 491.1;   Pc = 40.50000;  w = 0.3019;    %naphthalene
%Tc = 568.8;   Tb = 398.8;   Pc = 24.90000;  Zra = 0.2571;  %octane
%Tc = 692.502; Tb = 463.413; Pc = 31.60940;  w = 0.26945;   %decalin
%Tc = 658.2;   Tb = 489.5;   Pc = 18.20000;  Zra = 0.2466;  %dodecane
%Tc = 722;     Tb = 560;     Pc = 14.10000;  Zra = 0.2388;  %hexadecane

% Old code to be turned on/off for modification and testing purposes
%A = -7.63495; B = 1.50724; C = -3.19678; D = -2.78710; %p-xylene (1)
%A = -7.85178; B = 2.17172; C = -3.70505; D = -4.81238; %Naphthalene (1)
%s = 7.4936; beta = -4.41938; gamma = 0.17823; m = 1.47145; sigma0 = 55.72; mu = 1.30;
%A = 9.1832; B = 3631.55; C = -67.72; D = -1; %DECALIN RATIO GIVEN (3)
%s = 8.5739; beta = -4.31913; gamma = 0.29734; m = 1.63814; sigma0 = 54; mu = 1.28; 
%s = 9.7653; beta = -4.08617; gamma = 0.49323; m = 1.82630; sigma0 = 55.2; mu = 1.33; 





%%%%%%%%%%%% Calculation of Key Properties of each Component %%%%%%%%%%%%%




% There are six components in this surrogate mixture, and one of them,
% decalin, is split into its cis- and trans- forms. Key properties are
% calculated for all. 



% PARA-XYLENE 
% Key properties for para-xylene
Vstar_px = 0.3740;
wsrk_px = 0.3216;
sigma0 = 0;
mu = 0;
s = 0;
beta = 0;
gamma = 0;
m = 0;
A = 0;
B = 0;
C = 0;
D = 0;
Zra = 0;

% List key properties for para-xylene
Tc_px = 616.2;   Tb = 411.5;   Pc_px = 35.10000;  Zra = 0.2592;  
A = -7.63495; B = 1.50724; C = -3.19678; D = -2.78710; 

% Set para-xylene temperature equal to the same reduced temperature as the
% mixture.
T_px = Tc_px*(T_guess/Tcm); 

% Calculate surface tension for p-xylene
sigma_px = 63.21153*(1 - (T_px/Tc_px))^1.24262; % mN/m

% Calculate saturation pressure for para-xylene by calling subfunction
% pccalc
ps_px = pccalc(T_px,s,Tc_px,Pc_px,beta,gamma,m,A,B,C,D); % bar

% Old code for testing and modification
%vl_px = vlcalc(T_px,Vstar_px,Tc_px,wsrk_px); % m3/mol
%vl_px = vl_px/0.10616; % m3/kg





% NAPHTHALENE
% Key properties for naphthalene
Vstar_nap = 0.3834;
wsrk_nap = 0.3;
sigma0 = 0;
mu = 0;
s = 0;
beta = 0;
gamma = 0;
m = 0;
A = 0;
B = 0;
C = 0;
D = 0;
Zra = 0;

% More key properties for naphthalene
Tc_nap = 748.4;   Tb = 491.1;   Pc_nap = 40.50000;  w = 0.3019; 
A = -7.85178; B = 2.17172; C = -3.70505; D = -4.81238; 

% Same reduced temperature as the mixture
T_nap = Tc_nap*(T_guess/Tcm);

% Calculate surface tension for naphthalene
sigma_nap = 68.14*(1 - (T_nap/Tc_nap))^1.14593; % mN/m

% Calculate saturation pressure for naphthalene by calling subfunction
% pccalc
ps_nap = pccalc(T_nap,s,Tc_nap,Pc_nap,beta,gamma,m,A,B,C,D); % bar

%Old code for testing/modification
%vl_nap = vlcalc(T_nap,Vstar_nap,Tc_nap,wsrk_nap); % m3/mol
%vl_nap = vl_nap/0.1281705; % m3/kg





% OCTANE
% List key properties for Octane
Vstar_oct = 0.4904;
wsrk_oct = 0.3998;
sigma0 = 0;
mu = 0;
s = 0;
beta = 0;
gamma = 0;
m = 0;
A = 0;
B = 0;
C = 0;
D = 0;
Zra = 0;

% More key properties for octane
Tc_oct = 568.8;   Tb = 398.8;   Pc_oct = 24.90000;  Zra = 0.2571; 
s = 7.4936; beta = -4.41938; gamma = 0.17823; m = 1.47145; sigma0 = 55.72; mu = 1.30;

% Same reduced temperature as the mixture
T_oct = Tc_oct*(T_guess/Tcm);

% Calculate surface tension for octane
sigma_oct = 55.72*(1 - (T_oct/Tc_oct))^1.3; % mN/m

% Calculate saturation pressure for octane by calling subfunction pccalc
ps_oct = pccalc(T_oct,s,Tc_oct,Pc_oct,beta,gamma,m,A,B,C,D); % bar

% Old code for testing/modification
%vl_oct = vlcalc(T_oct,Vstar_oct,Tc_oct,wsrk_oct); % m3/mol
%vl_oct = vl_oct/0.11423; % m3/kg





% CIS - DECALIN
% Key properties for cis-decalin component
sigma0 = 0;
mu = 0;
s = 0;
beta = 0;
gamma = 0;
m = 0;
A = 0;
B = 0;
C = 0;
D = 0;
Zra = 0;

% More key properties for cis-decalin
Tc_cdec = 702.3; Tb = 468.9; Pc_cdec = 32;  w = 0.286; 
A = 9.211; B = 3671.61; C = -69.74; D = -1; 

% Calculate acentric factor - an additional factor for cis-decalin (needs
% to be done for trans-decalin as well). 
wsrk_cdec = 0.286;
R = 0.083144;
Vstar_cdec = (R*Tc_cdec/Pc_cdec)*(0.2717636 - 0.05759377*wsrk_cdec + ...
    0.05527757*(wsrk_cdec^2));

% Same reduced temperature as the mixture
T_cdec = Tc_cdec*(T_guess/Tcm);

% Calculate surface tension for cis-decalin
sigma_cdec = 64.62*(1 - (T_cdec/Tc_cdec))^1.2925; % mN/m

% Calculate saturation pressure for cis-decalin
x_cdec = 1 - (T_cdec/Tc_cdec);
ps_cdec = Pc_cdec*(1 - 7.1508*x_cdec + 19.59405*x_cdec^2 ...
    - 24.60417*x_cdec^3 + 12.068*x_cdec^4); % bar

% Old code for testing/modification
%vl_cdec = vlcalc(T_cdec,Vstar_cdec,Tc_cdec,wsrk_cdec); % m3/mol
%vl_cdec = vl_cdec/0.13825; % m3/kg





% TRANS - DECALIN
% List key properties for trans-decalin
sigma0 = 0;
mu = 0;
s = 0;
beta = 0;
gamma = 0;
m = 0;
A = 0;
B = 0;
C = 0;
D = 0;
Zra = 0;

% More key properties for trans-decalin
Tc_tdec = 687.1; Tb = 460.5; Pc_tdec = 31.4;  w = 0.270;
A = 9.1787; B = 3610.66; C = -66.49; D = -1;

% acentric factor prediction
wsrk_tdec = 0.270;
R = 0.083144;
Vstar_tdec = (R*Tc_tdec/Pc_tdec)*(0.2717636 - 0.05759377*wsrk_tdec + ...
    0.05527757*(wsrk_tdec^2));

% Same reduced temperature as the mixture
T_tdec = Tc_tdec*(T_guess/Tcm);

% Calculate surface tension for trans-decalin
sigma_tdec = 61.34*(1 - (T_tdec/Tc_tdec))^1.29579; % mN/m

% Calculate saturation pressure for trans-decalin
x_tdec = 1 - (T_tdec/Tc_tdec);
ps_tdec = Pc_tdec*(1 - 6.64617*x_tdec + 15.6026*x_tdec^2 ...
    - 14.19751*x_tdec^3 + 3.11617*x_tdec^4); % bar

% Old code for testing/modification
%vl_tdec = vlcalc(T_tdec,Vstar_tdec,Tc_tdec,wsrk_tdec); % m3/mol
%vl_tdec = vl_tdec/0.13825; % m3/kg





% DODECANE
% List key properties for dodecane
Vstar_dod = 0.7558;
wsrk_dod = 0.5807;
sigma0 = 0;
mu = 0;
s = 0;
beta = 0;
gamma = 0;
m = 0;
A = 0;
B = 0;
C = 0;
D = 0;
Zra = 0;

% List more properties for dodecane
Tc_dod = 658.2;   Tb = 489.5;   Pc_dod = 18.20000;  Zra = 0.2466;
s = 8.5739; beta = -4.31913; gamma = 0.29734; m = 1.63814; sigma0 = 54; mu = 1.28;

% Same reduced temperature as the mixture
T_dod = Tc_dod*(T_guess/Tcm);

% Calculate surface tension for dodecane
sigma_dod = 54*(1 - (T_dod/Tc_dod))^1.28; % mN/m

% Calculate saturation pressure for dodecane by calling subfunction pccalc
ps_dod = pccalc(T_dod,s,Tc_dod,Pc_dod,beta,gamma,m,A,B,C,D); % bar

% Old code for testing/modification
%vl_dod = vlcalc(T_dod,Vstar_dod,Tc_dod,wsrk_dod); % m3/mol
%vl_dod = vl_dod/0.17034; % m3/kg





% HEXADECANE
% List key properties for hexadecane
Vstar_hex = 1.0539;
wsrk_hex = 0.7667;
sigma0 = 0;
mu = 0;
s = 0;
beta = 0;
gamma = 0;
m = 0;
A = 0;
B = 0;
C = 0;
D = 0;
Zra = 0;

% More properties for hexadecane
Tc_hex = 722;     Tb = 560;     Pc_hex = 14.10000;  Zra = 0.2388; 
s = 9.7653; beta = -4.08617; gamma = 0.49323; m = 1.82630; sigma0 = 55.2; mu = 1.33;

% Same reduced temperature as the mixture
T_hex = Tc_hex*(T_guess/Tcm);

% Calculate surface tension of hexadecane
sigma_hex = 55.2*(1 - (T_hex/Tc_hex))^1.33; % mN/m

% Calculate saturation pressure for hexadecane by calling subfunction
% pccalc
ps_hex = pccalc(T_hex,s,Tc_hex,Pc_hex,beta,gamma,m,A,B,C,D); % bar

% Old code for testing/modification
%vl_hex = vlcalc(T_hex,Vstar_hex,Tc_hex,wsrk_hex); % m3/mol
%vl_hex = vl_hex/0.226448; % m3/kg





%%%%%%%%%%%%%% Code for Calculation of Mixture Properties %%%%%%%%%%%%%%%%




% Get mole fractions for the mixture
% vol frac * density / molecular weight = # mols of comp
% https://www.atsdr.cdc.gov/toxprofiles/tp67-c4.pdf

% Calculation of component molar amounts
m1_px = 0.085*861/0.10616; 
m1_nap = 0.08*1145/0.1281705; 
m1_oct = 0.035*703/0.11423;
m1_cdec = 0.35*0.4*896/0.13825;
m1_tdec = 0.35*0.6*896/0.13825;
m1_dod = 0.4*750/0.17034;
m1_hex = 0.05*770/0.22645;

% Total number of moles of the surrogate mixture
m_total = m1_px + m1_nap + m1_oct + m1_cdec + m1_tdec + m1_dod + m1_hex; 

% Calculate component mole fractions
mf_px = m1_px/m_total;
mf_nap = m1_nap/m_total;
mf_oct = m1_oct/m_total;
mf_cdec = m1_cdec/m_total;
mf_tdec = m1_tdec/m_total;
mf_dod = m1_dod/m_total;
mf_hex = m1_hex/m_total;


%Calculate surface tension of the mixture
sigma_m = mf_px*sigma_px + mf_nap*sigma_nap + mf_oct*sigma_oct...
    + mf_cdec*sigma_cdec + mf_tdec*sigma_tdec + mf_dod*sigma_dod + mf_hex*sigma_hex;
sigma_m = sigma_m*0.001; % mN/m to N/m


%Calculate saturation pressure of the mixture
ps_m = (mf_px*ps_px/Pc_px + mf_nap*ps_nap/Pc_nap + mf_oct*ps_oct/Pc_oct...
    + mf_cdec*ps_cdec/Pc_cdec + mf_tdec*ps_tdec/Pc_tdec...
    + mf_dod*ps_dod/Pc_dod + mf_hex*ps_hex/Pc_hex)*Pcm; %bar


% Old code for testing and validation purposes
%vl_m = mf_px*vl_px + mf_nap*vl_nap + mf_oct*vl_oct + mf_cdec*vl_cdec +...
%mf_tdec*vl_tdec + mf_dod*vl_dod + mf_hex*vl_hex; %m3/kg


% Calculation of interaction pair constants

% para-xylene and the six
comp_px_px = mf_px*mf_px*(Vstar_px*Tc_px*Vstar_px*Tc_px)^(1/2);
comp_px_nap = mf_px*mf_nap*(Vstar_px*Tc_px*Vstar_nap*Tc_nap)^(1/2);
comp_px_oct = mf_px*mf_oct*(Vstar_px*Tc_px*Vstar_oct*Tc_oct)^(1/2);
comp_px_cdec = mf_px*mf_cdec*(Vstar_px*Tc_px*Vstar_cdec*Tc_cdec)^(1/2);
comp_px_tdec = mf_px*mf_tdec*(Vstar_px*Tc_px*Vstar_tdec*Tc_tdec)^(1/2);
comp_px_dod = mf_px*mf_dod*(Vstar_px*Tc_px*Vstar_dod*Tc_dod)^(1/2);
comp_px_hex = mf_px*mf_hex*(Vstar_px*Tc_px*Vstar_hex*Tc_hex)^(1/2);

% naphthalene and the six
comp_nap_px = mf_nap*mf_px*(Vstar_nap*Tc_nap*Vstar_px*Tc_px)^(1/2);
comp_nap_nap = mf_nap*mf_nap*(Vstar_nap*Tc_nap*Vstar_nap*Tc_nap)^(1/2);
comp_nap_oct = mf_nap*mf_oct*(Vstar_nap*Tc_nap*Vstar_oct*Tc_oct)^(1/2);
comp_nap_cdec = mf_nap*mf_cdec*(Vstar_nap*Tc_nap*Vstar_cdec*Tc_cdec)^(1/2);
comp_nap_tdec = mf_nap*mf_tdec*(Vstar_nap*Tc_nap*Vstar_tdec*Tc_tdec)^(1/2);
comp_nap_dod = mf_nap*mf_dod*(Vstar_nap*Tc_nap*Vstar_dod*Tc_dod)^(1/2);
comp_nap_hex = mf_nap*mf_hex*(Vstar_nap*Tc_nap*Vstar_hex*Tc_hex)^(1/2);

% octane and the six
comp_oct_px = mf_oct*mf_px*(Vstar_oct*Tc_oct*Vstar_px*Tc_px)^(1/2);
comp_oct_nap = mf_oct*mf_nap*(Vstar_oct*Tc_oct*Vstar_nap*Tc_nap)^(1/2);
comp_oct_oct = mf_oct*mf_oct*(Vstar_oct*Tc_oct*Vstar_oct*Tc_oct)^(1/2);
comp_oct_cdec = mf_oct*mf_cdec*(Vstar_oct*Tc_oct*Vstar_cdec*Tc_cdec)^(1/2);
comp_oct_tdec = mf_oct*mf_tdec*(Vstar_oct*Tc_oct*Vstar_tdec*Tc_tdec)^(1/2);
comp_oct_dod = mf_oct*mf_dod*(Vstar_oct*Tc_oct*Vstar_dod*Tc_dod)^(1/2);
comp_oct_hex = mf_oct*mf_hex*(Vstar_oct*Tc_oct*Vstar_hex*Tc_hex)^(1/2);

% cis-decalin and the six
comp_cdec_px = mf_cdec*mf_px*(Vstar_cdec*Tc_cdec*Vstar_px*Tc_px)^(1/2);
comp_cdec_nap = mf_cdec*mf_nap*(Vstar_cdec*Tc_cdec*Vstar_nap*Tc_nap)^(1/2);
comp_cdec_oct = mf_cdec*mf_oct*(Vstar_cdec*Tc_cdec*Vstar_oct*Tc_oct)^(1/2);
comp_cdec_cdec = mf_cdec*mf_cdec*(Vstar_cdec*Tc_cdec*Vstar_cdec*Tc_cdec)^(1/2);
comp_cdec_tdec = mf_cdec*mf_tdec*(Vstar_cdec*Tc_cdec*Vstar_tdec*Tc_tdec)^(1/2);
comp_cdec_dod = mf_cdec*mf_dod*(Vstar_cdec*Tc_cdec*Vstar_dod*Tc_dod)^(1/2);
comp_cdec_hex = mf_cdec*mf_hex*(Vstar_cdec*Tc_cdec*Vstar_hex*Tc_hex)^(1/2);

% trans-decalin and the six
comp_tdec_px = mf_tdec*mf_px*(Vstar_tdec*Tc_tdec*Vstar_px*Tc_px)^(1/2);
comp_tdec_nap = mf_tdec*mf_nap*(Vstar_tdec*Tc_tdec*Vstar_nap*Tc_nap)^(1/2);
comp_tdec_oct = mf_tdec*mf_oct*(Vstar_tdec*Tc_tdec*Vstar_oct*Tc_oct)^(1/2);
comp_tdec_cdec = mf_tdec*mf_cdec*(Vstar_tdec*Tc_tdec*Vstar_cdec*Tc_cdec)^(1/2);
comp_tdec_tdec = mf_tdec*mf_tdec*(Vstar_tdec*Tc_tdec*Vstar_tdec*Tc_tdec)^(1/2);
comp_tdec_dod = mf_tdec*mf_dod*(Vstar_tdec*Tc_tdec*Vstar_dod*Tc_dod)^(1/2);
comp_tdec_hex = mf_tdec*mf_hex*(Vstar_tdec*Tc_tdec*Vstar_hex*Tc_hex)^(1/2);

% dodecane and the six
comp_dod_px = mf_dod*mf_px*(Vstar_dod*Tc_dod*Vstar_px*Tc_px)^(1/2);
comp_dod_nap = mf_dod*mf_nap*(Vstar_dod*Tc_dod*Vstar_nap*Tc_nap)^(1/2);
comp_dod_oct = mf_dod*mf_oct*(Vstar_dod*Tc_dod*Vstar_oct*Tc_oct)^(1/2);
comp_dod_cdec = mf_dod*mf_cdec*(Vstar_dod*Tc_dod*Vstar_cdec*Tc_cdec)^(1/2);
comp_dod_tdec = mf_dod*mf_tdec*(Vstar_dod*Tc_dod*Vstar_tdec*Tc_tdec)^(1/2);
comp_dod_dod = mf_dod*mf_dod*(Vstar_dod*Tc_dod*Vstar_dod*Tc_dod)^(1/2);
comp_dod_hex = mf_dod*mf_hex*(Vstar_dod*Tc_dod*Vstar_hex*Tc_hex)^(1/2);

% hexadecane and the six
comp_hex_px = mf_hex*mf_px*(Vstar_hex*Tc_hex*Vstar_px*Tc_px)^(1/2);
comp_hex_nap = mf_hex*mf_nap*(Vstar_hex*Tc_hex*Vstar_nap*Tc_nap)^(1/2);
comp_hex_oct = mf_hex*mf_oct*(Vstar_hex*Tc_hex*Vstar_oct*Tc_oct)^(1/2);
comp_hex_cdec = mf_hex*mf_cdec*(Vstar_hex*Tc_hex*Vstar_cdec*Tc_cdec)^(1/2);
comp_hex_tdec = mf_hex*mf_tdec*(Vstar_hex*Tc_hex*Vstar_tdec*Tc_tdec)^(1/2);
comp_hex_dod = mf_hex*mf_dod*(Vstar_hex*Tc_hex*Vstar_dod*Tc_dod)^(1/2);
comp_hex_hex = mf_hex*mf_hex*(Vstar_hex*Tc_hex*Vstar_hex*Tc_hex)^(1/2);


% Calculation of key molar physical property of the mixture.
Vstartot1 = mf_px*Vstar_px + mf_nap*Vstar_nap + mf_oct*Vstar_oct...
    + mf_cdec*Vstar_cdec + mf_tdec*Vstar_tdec + mf_dod*Vstar_dod...
    + mf_hex*Vstar_hex;

Vstartot2 = mf_px*(Vstar_px^(2/3)) + mf_nap*(Vstar_nap^(2/3))...
    + mf_oct*(Vstar_oct^(2/3)) + mf_cdec*(Vstar_cdec^(2/3))...
    + mf_tdec*(Vstar_tdec^(2/3)) + mf_dod*(Vstar_dod^(2/3))...
    + mf_hex*(Vstar_hex^(2/3));

Vstartot3 = mf_px*(Vstar_px^(1/3)) + mf_nap*(Vstar_nap^(1/3))...
    + mf_oct*(Vstar_oct^(1/3)) + mf_cdec*(Vstar_cdec^(1/3))...
    + mf_tdec*(Vstar_tdec^(1/3)) + mf_dod*(Vstar_dod^(1/3))...
    + mf_hex*(Vstar_hex^(1/3));

% Calculation of the key property: Mixture V*
Vstar_m = (1/4)*(Vstartot1 + 3*Vstartot2*Vstartot3);

%Calculation of the liquid-molar critical temperature of the mixture
Tcvl_m = (comp_px_px + comp_px_nap + comp_px_oct + comp_px_cdec...
    + comp_px_tdec + comp_px_dod + comp_px_hex + comp_nap_px...
    + comp_nap_nap + comp_nap_oct + comp_nap_cdec + comp_nap_tdec...
    + comp_nap_dod + comp_nap_hex + comp_oct_px + comp_oct_nap...
    + comp_oct_oct + comp_oct_cdec + comp_oct_tdec + comp_oct_dod...
    + comp_oct_hex + comp_cdec_px + comp_cdec_nap + comp_cdec_oct...
    + comp_cdec_cdec + comp_cdec_tdec + comp_cdec_dod + comp_cdec_hex...
    + comp_tdec_px + comp_tdec_nap + comp_tdec_oct + comp_tdec_cdec...
    + comp_tdec_tdec + comp_tdec_dod + comp_tdec_hex + comp_dod_px...
    + comp_dod_nap + comp_dod_oct + comp_dod_cdec + comp_dod_tdec...
    + comp_dod_dod + comp_dod_hex + comp_hex_px + comp_hex_nap...
    + comp_hex_oct + comp_hex_cdec + comp_hex_tdec + comp_hex_dod...
    + comp_hex_hex)/Vstar_m;

% Old code for verification purposes and testing
%Tcvl_m = Tcm;



% Calculation of the acentric factor for the mixture - key physical prop.
wsrk_m = mf_px*wsrk_px + mf_nap*wsrk_nap + mf_oct*wsrk_oct...
    + mf_cdec*wsrk_cdec + mf_tdec*wsrk_tdec + mf_dod*wsrk_dod...
    + mf_hex*wsrk_hex;

% Calculation of liquid-molar volume of the mixture by calling subfunction
% vlcalc
vl_m = vlcalc(T_guess,Vstar_m,Tcvl_m,wsrk_m); %m3/mol

% Mixture molecular mass
m_mol_mass = mf_px*0.10616 + mf_nap*0.1281705 + mf_oct*0.11423...
    + mf_cdec*0.13825 + mf_tdec*0.13825 + mf_dod*0.17034 + mf_hex*0.22645;

% Conversion of liquid molar volume to mass-based form
vl_m = vl_m/m_mol_mass; %m3/kg









%%%%%%%%%%%%%%%%%%%%%%%%% Subfunction sigmacalc %%%%%%%%%%%%%%%%%%%%%%%%%%
    function sigma = sigmacalc(T,sigma0,mu,Tc,Tb,Pc)
        
        % Reduced temperature and boiling temperature
        Tr = T/Tc;
        Tbr = Tb/Tc;
        
        % Sort through component properties with "if" statement test
        if sigma0 ~= 0
            sigma = sigma0*(1 - Tr)^mu;
        else
            % For napthalane, p-xylene, decalin
            Q = 0.1196*(1 + ((Tbr*log(Pc/1.01325))/(1 - Tbr))) - 0.279;
            sigma = ((Pc)^(2/3))*((Tc)^(1/3))*Q*((1-Tr)^(11/9));
            % Surface tension (sigma) (from Properties of Gases & Liquids,
            % Reid 1987, p. 637) (Brock and Bird Equation)
        end
    end
%%%%%%%%%%%%%%%%%%%%% End of Subfunction sigmacalc %%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% Subfunction pccalc %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Ps = pccalc(T,s,Tc,Pc,beta,gamma,m,A,B,C,D)
        
        % Reduced temperature of the chemical
        Tr = T/Tc;
        
        % Chemical property
        n = 7;
        
        % Sift through chemical properties in "if" statement to match the
        % correct form of the equation to the proper chemicals
        % respectively.
        if s ~= 0
            
            % Calculate saturation pressure
            Pc = Pc/1.01325; %convert input Pc from bar to atm
            Ps = Pc*exp(beta*((1/((Tr)^m))-1) + gamma*(((Tr)^n) - 1));
            Ps = Ps*1.01325; %atm to bar
        
        elseif D == -1
            % Calculate saturation pressure
            Ps = exp(A - (B/(T + C)));
        
        else
            % Calculate saturation pressure
            x1 = 1 - Tr;
            Ps = Pc*exp(((1-x1)^(-1))*((A*x1) + B*(x1^1.5) + C*(x1^3) + D*(x1^6))); %bar
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%% End of Subfunction pccalc %%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% Subfunction vlcalc %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function vl = vlcalc(T,Vstar,Tc,wsrk)
        
        % Definition of key constants for formulas
        R = 0.083144;
        a = -1.52816;    b = 1.43907;    c = -0.81446;    d = 0.190454;
        e = -0.296123;    f = 0.386914;    g = -0.0427258;    h = -0.0480645;
        
        % Reduced temperature
        Tr = T/Tc;
        
        % Calculation of key constants for equation
        VR0 = 1 + (a*(1-Tr)^(1/3)) + (b*(1-Tr)^(2/3)) + (c*(1-Tr))...
            + (d*(1-Tr)^(4/3));
        VRd = (e + f*Tr + g*(Tr^2) + h*(Tr^3))/(Tr - 1.00001);
        
        % Computation equation for liquid-molar volume of chemical
        vl = Vstar*VR0*(1 - wsrk*VRd);
        vl = vl/1000; %m3/mol    
    end
%%%%%%%%%%%%%%%%%%%%%% End of Subfunction vlcalc %%%%%%%%%%%%%%%%%%%%%%%%%


end

