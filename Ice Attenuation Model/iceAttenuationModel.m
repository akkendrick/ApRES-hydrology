function [ T_cel NgreenTot ] = iceAttenuationModel()
% Compute ice attenuation model based on series of papers by MacGregor et
% al and Matsouka et al. 
% Outputs:
%       Currently set to give the one way attenuation rate of Greenland ice
%       in dB/km from -50 to 0 degrees C
% Alex Kendrick
% 1/30/2017


%close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute radar power

% From Matsuoka et al. 2012
% Pdb = Sdb + Rdb - Ldb - Bdb

%Uniform permittivity from Matsuoka (2010)
epsilon = 3.2; % Equivalent to 168 m/us propagation speed
f = 300; %Center frequency is 600 MHz

% From Matsuoka et al. 2010 (2) 
% i = 1: Water i = 2: Acidity i = 3: Salinity
% Parameters from MacGregor et al 2007 Table 2
% Molar conductivity of H+ : uHcl = 3.2 S/(m*M)
% Molar conductivity of Cl- : uCl = 0.43 S/(m*M)
% Pure ice conductivity: 6.6 uS/m

% Activation Energies
% Epure = 0.55 eV
% EH+ = 0.20 eV
% ECl- = 0.19 eV

Tr = 251; %Ref Temp in K
T_cel = -35:0.001:15;
% Greenland bed depth (taken from MIMO processing)
H = 617.23;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Greenland Params from Macgregor 2015 M07 model, Holocene epoch
Tr_green = -21+273;
C_green = [1 1.6 0.4 0.5];
sigma_green = [9.2 3.2 0.43 0.8];
E_green = [0.51 0.20 0.19 0.23];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate attenuation
%C = [1 1.3 4.2]; %Matsouka 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matsouka 2012
C_siple = [ 1 1.2 4.1]; 
sigma_siple = [9.2 3.2 0.43]; 
E_siple = [0.51 0.20 0.19];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%C = [1 1.3 4.2]; % Concentrations in uM (mol/uL)
%sigma = [6.6 3.2 0.43]; % Conductivities in 
%E = [0.55 0.20 0.19]; % Activation energies in eV
T = T_cel + 273; % Ice temperature *** This is a function of depth
k = 8.61733*10^-5; % Boltzmann const in ev/K

for j = 1:length(T)
    
    for i = 1:3
        Nsiple(i) = sigma_siple(i)*C_siple(i)*exp(-E_siple(i)/k*(1/T(j) - 1/Tr));
    end
    
    for i =1:4
       Ngreen(i) = sigma_green(i)*C_green(i)*exp(-E_green(i)/k*(1/T(j) - 1/Tr));
    end
        
    NgreenTot(j) = 0.914*sum(Ngreen);
    
    NsipleTot(j) = 0.914*sum(Nsiple);
    %NdB(j) = 10*log10(1000*(j));

    Npure(j) = 0.914*sigma_siple(1)*C_siple(1)*exp(-(E_siple(1)/k) * (1/T(j) - 1/Tr));
    %Npure(j) = 10*log10(1000*Npure(j));

    NGud(j) = 0.914*(15.4*exp(-0.33/k*(1/T(j) - 1/Tr)));
    %NGud(j) = 10*log10(1000*NGud(j));
    
end

figure(10)
hold on
plot(T_cel, NsipleTot,'LineWidth',2)
plot(T_cel, NGud,'LineWidth',2)
plot(T_cel, Npure,'LineWidth',2)
plot(T_cel, NgreenTot,'LineWidth',2)

xlabel('Temperature (\circC)')
ylabel('Attenuation rate N (dB/km, one way')

legend('Siple Dome','Gudmandsen','Pure Ice','Greenland')
% Calculate Fresnel Reflectivity
%R = 20*log10(4.5*10^-3/f*

% h = 1;
% Gdb = 2*10*log10(h + z2/sqrt(epsilon))


end

