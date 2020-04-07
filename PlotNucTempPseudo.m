function [PseudoNuc] = PlotNucTempPseudo()
% This is the main host function of this jet fuel surrogate program bundle.
% The goal of this program bundle is to predict the thermodynamic
% conditions under which a given synthetic jet fuel mixture, or surrogate,
% releases the most energy upon combustion. A goal is to calculate and plot
% nucleation temperatures of the chosen surrogate mixuture in order to
% predict its superheat limit, a fundamental temperature that is linked to
% optimal fuel performance. 

% This function plots the nucleation temperatures for our selected
% six-component surrogate mixture against the nucleation rate. It treats 
% the mixture as if it were one component, hence, a
% "pseudo-single-component mixture". The six components are para-xylene,
% naphthalene, octane, decalin, dodecane, and hexadecane.

% With this nucleation temperature curve, assymptotes can be constructed
% either mathematically in the form of a limit, or graphically using a
% trial-and-error-based method of prediction. These asymptotes can be used
% to estimate the superheat limit of the jet fuel surrogate.





%%%%%%%%%%%%%%%%%%%%%%%%% Set Up Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%




%guess Tnuc estimate to run initial runthrough of temperature computation.
T_guess = 600; 

% List important physical constants.
K = 1.38064852e-23; %Boltsmann Constant, ESI units
Av = 6.0221409e23; %Avagadro's Number
Po = 1.01325; %atmospheric pressure in bar
R = 8.314e-5; %gas constant, m3-bar/mole-K
J = 10^21; %"Guess" J value

% Calculate molar amounts of the six components (cis- and trans- decalin
% together form the decalin component). 
m1_px = 0.085*861/0.10616;
m1_nap = 0.08*1145/0.1281705;
m1_oct = 0.035*703/0.11423;
m1_cdec = 0.35*0.4*896/0.13825;
m1_tdec = 0.35*0.6*896/0.13825;
m1_dod = 0.4*750/0.17034;
m1_hex = 0.05*770/0.22645;

% Calculate total number of moles of the mixture 
m_total = m1_px + m1_nap + m1_oct + m1_cdec + m1_tdec + m1_dod + m1_hex;

% Calculate mole fractions of the different components
mf_px = m1_px/m_total;
mf_nap = m1_nap/m_total;
mf_oct = m1_oct/m_total;
mf_cdec = m1_cdec/m_total;
mf_tdec = m1_tdec/m_total;
mf_dod = m1_dod/m_total;
mf_hex = m1_hex/m_total;

% Calculate molecular weight fractions for each of the components.
wf_px = mf_px*0.10616;
wf_nap = mf_nap*0.1281705;
wf_oct = mf_oct*0.11423;
wf_cdec = mf_cdec*0.13825;
wf_tdec = mf_tdec*0.13825;
wf_dod = mf_dod*0.17034;
wf_hex = mf_hex*0.22645;

% Calculate total molecular weight of the mixture
W = wf_px + wf_nap + wf_oct + wf_cdec + wf_tdec + wf_dod + wf_hex;





%%%%%% Iterate until convergence to plot nucleation temperatures %%%%%%%%




% Cycle through sample nucleation rates J. These will be the points where
% nucleation temperatures for the surrogate mixture are plotted against
% nucleation rate.
for Jcount = 1:1:21
    J(Jcount) = (10^21)*(10^(1-Jcount));
    T_guess = 627;
    
    % iterate many times beginning with the initial "guessed" nucleation
    % temperature of the surrogate, until convergence is reached between
    % each new guess and the nucleation temperature calculated based on
    % that "guess" as an input.
    for counter = 1:2000
        T_guess = T_guess + 0.01;
        NucTemp(counter) = T_guess;
        
        % three critical properties (mixturesurface tension sigma_m,
        % mixture saturation pressure ps_m, and mixture liquid-molar volume
        % vl_m) are calculated for the mixture based on the guessed
        % nulceation temperature T_guess. This is the temperature that is
        % iterating through until convergence with the calculated surrogate
        % nucleation temperature below. This calls another function in this
        % bundle to do so.
        [sigma_m,ps_m,vl_m] = get3propertiesnew3(T_guess);
        P = ps_m*exp((vl_m*W/(R*T_guess))*((Po - ps_m)^1));
        delta_omega = (10e-10)*16*pi*(sigma_m^3)/(3*((P-Po)^2));
        
        % Calculation of constants for nucleation temperature calculation.
        m = W/Av;
        
        No = (Av/(W*vl_m))^(2/3);
        
        C = No*sqrt(2*sigma_m/(pi*m));
        
        % Difference (to be reduced to convergence) between iterated
        % guessed nucleation temperature, and the calculated surogate
        % nucleation temperature below.
        Cdiff(counter) = T_guess - delta_omega*((log(C/J(Jcount)))^-1)/K;
    end
    
    % Clear convergence figure
    figure(2); clf;
    
    [~,idx] = min(abs(Cdiff));
    
    Nuc_Temp(Jcount) = NucTemp(idx);
    J(Jcount) = log10(J(Jcount));
    
end






%%%%%%%%%%%%%%%%% Plot and post process nucleation data %%%%%%%%%%%%%%%%%%




% plot surrogate nucleation temperatures against the selected surrogate
% nucleation rates.
figure(2); clf;
plot(J,Nuc_Temp,'k*');

% Old Code to be turned on/off
%axis([5,9,520,630]);

% Label Axes of plot
xlabel('Logarithmic Surrogate Nucleation Rate (nuclei/m^2-s)');
ylabel('Surrogate Superheat Limit (K)');
title('Surrogate Superheat Limit');

% Old Code to be turned on/off
%blank = zeros(length(Cdiff));
%plot(blank,'b--');


%Format output data for surrogate nucleation temperature and surrogate
%nucleation rate.
J = J';
Nuc_Temp = Nuc_Temp';
Jnew = zeros(1,21);
Nuc_Tempnew = zeros(1,21);

for counter1 = 1:21
Jnew(22-counter1) = J(counter1);
Nuc_Tempnew(22-counter1) = Nuc_Temp(counter1);
end

% matrix that stores surrogate nucleation temperature and surrogate
%nucleation rate.
PseudoNuc(:,1) = Jnew';
PseudoNuc(:,2) = Nuc_Tempnew';


end

