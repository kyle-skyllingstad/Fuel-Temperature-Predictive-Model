function [Tcm,Pcm] = MixtureCriticalTP()
% This function calculates the critical temperature and pressure for the
% surrogate mixture to be used in further calculations. This function is
% called by the other two functions in this bundle.





%%%%%%%%%%%%%%%%%%%%%%%%% Set Up Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Calculation of constants for formulas
R = 8.314e-5; %gas constant, m3-bar/mole-K

% Definition of critical pressure for the six components (decalin is split
% into its cis- and trans- forms.
Pc_px = 35.10000;   %p-xylene
Pc_nap = 40.50000;  %naphthalene
Pc_oct = 24.90000;  %octane
Pc_cdec = 32;  %cis-decalin
Pc_tdec = 31.4;  %trans-decalin
Pc_dod = 18.20000;  %dodecane
Pc_hex = 14.10000;  %hexadecane

% Definition of the critical temperatures of the six components (decalin
% split). 
Tc_px = 616.2;
Tc_nap = 748.4;
Tc_oct = 568.8;
Tc_cdec = 702.3;
Tc_tdec = 687.1;
Tc_dod = 658.2;
Tc_hex = 722;

% Literature-based physical constants calculated for the six components, to
% be used in further calculation.
Zc_px = 0.26;
Zc_nap = 0.269;
Zc_oct = 0.259;
Zc_cdec = 0.250; 
Zc_tdec = 0.254; 
Zc_dod = 0.24;
Zc_hex = 0.2142;

%Critical volumes calculated for the six components.
Vc_px = 379;
Vc_nap = 413;
Vc_oct = 492;
Vc_cdec = 530;  % 0.4 ratio cis-
Vc_tdec = 530;  % 0.6 ratio trans-
Vc_dod = 713;
Vc_hex = 678;    

% Critical references from which these data were obtained.
% for Zc: www.egr.msu.edu/~lira/computer/EXCEL/PREOS.XLS
% for more, go to http://www.eng.umd.edu/~nsw/chbe250/critical.dat

% Volume fractions for the six components.
volf_px = 0.085;
volf_nap = 0.08;
volf_oct = 0.035;
volf_cdec = 0.35*0.4;
volf_tdec = 0.35*0.6;
volf_dod = 0.4;
volf_hex = 0.05;

% calculate molar amounts of the six components
m1_px = 0.085*861/0.10616;
m1_nap = 0.08*1145/0.1281705;
m1_oct = 0.035*703/0.11423;
m1_cdec = 0.35*0.4*896/0.13825;
m1_tdec = 0.35*0.6*896/0.13825;
m1_dod = 0.4*750/0.17034;
m1_hex = 0.05*770/0.22645;

% Calculate the total number of moles of the mixture
m_total = m1_px + m1_nap + m1_oct + m1_cdec + m1_tdec + m1_dod + m1_hex;

% calculate mole fractions of each ofthe six components of the mixture.
mf_px = m1_px/m_total;
mf_nap = m1_nap/m_total;
mf_oct = m1_oct/m_total;
mf_cdec = m1_cdec/m_total;
mf_tdec = m1_tdec/m_total;
mf_dod = m1_dod/m_total;
mf_hex = m1_hex/m_total;

% Calculate modified volume total of the mixture
phitot = mf_px*Vc_px + mf_nap*Vc_nap + mf_oct*Vc_oct + mf_cdec*Vc_cdec...
    + mf_tdec*Vc_tdec + mf_dod*Vc_dod + mf_hex*Vc_hex;

% Calculate modified volume fractions of the six components
phi_px = mf_px*Vc_px/phitot;
phi_nap = mf_nap*Vc_nap/phitot;
phi_oct = mf_oct*Vc_oct/phitot;
phi_cdec = mf_cdec*Vc_cdec/phitot;
phi_tdec = mf_tdec*Vc_tdec/phitot;
phi_dod = mf_dod*Vc_dod/phitot;
phi_hex = mf_hex*Vc_hex/phitot;





%%%%%%%%%%%%%%%%%%%%%%%%% Calculation Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Calculate actual critical temperature of the mixture
Tcmactual = phi_px*Tc_px + phi_nap*Tc_nap + phi_oct*Tc_oct...
    + phi_cdec*Tc_cdec + phi_tdec*Tc_tdec + phi_dod*Tc_dod + phi_hex*Tc_hex;

% Calculate pseudo and usable critical temperature of the mixture
Tcm = mf_px*Tc_px + mf_nap*Tc_nap + mf_oct*Tc_oct + mf_cdec*Tc_cdec...
    + mf_tdec*Tc_tdec + mf_dod*Tc_dod + mf_hex*Tc_hex;




% Old code to be turned on/off for mixture critical pressure calculation
% Pcm = mf_px*Pc_px + mf_nap*Pc_nap + mf_oct*Pc_oct...
%     + mf_dec*Pc_dec + mf_dod*Pc_dod + mf_hex*Pc_hex;
% Pcm = 0.085*Pc_px + 0.08*Pc_nap + 0.035*Pc_oct + 0.35*Pc_dec + ...
% 0.4*Pc_dod + 0.05*Pc_hex; %bar


% Calculation of the mixture critical pressure
Pcm = 1000000*R*Tcm*(mf_px*Zc_px + mf_nap*Zc_nap + mf_oct*Zc_oct...
    + mf_cdec*Zc_cdec + mf_tdec*Zc_tdec + mf_dod*Zc_dod + mf_hex*Zc_hex)...
    /(mf_px*Vc_px + mf_nap*Vc_nap + mf_oct*Vc_oct + mf_cdec*Vc_cdec...
    + mf_tdec*Vc_tdec + mf_dod*Vc_dod + mf_hex*Vc_hex); %bar


% Old code to be turned on/off for mixture critical pressure calculation
%Pcm = mf_px*Pc_px + mf_nap*Pc_nap + mf_oct*Pc_oct + mf_cdec*Pc_cdec...
%+ mf_tdec*Pc_tdec + mf_dod*Pc_dod + mf_hex*Pc_hex

end

