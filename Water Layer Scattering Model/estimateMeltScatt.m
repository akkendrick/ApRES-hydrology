function [ T_water, T_water_iceMix ] = estimateMeltScatt(deltaAtten,scattFrac,fc, sigma_w,phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimateMeltScatt
% For a given observed change in attenuation estimate the amount of melt
% water needed to produce this change
%
% Inputs:
%   deltaAtten - change in attenuation
%   scattFrac - proportion of the attenuation attributable to volume
%       scattering
%   fc - center frequency
%   sigma_w - conductivity of water
%   phi - porosity
%
% Outputs:
%   T_water - estimated thickness of stored water
%   T_water_iceMix - estimated thickness of macroporous storage
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_o = 1.2566370614*10^(-6); %H/m;

% Conductivity mixing model and mean value of s from Geldsetzer et al 2009
s = 1.67;
%s = 1.00;

% The water ice mixture (phi or v) is considered to be the thickness of the 
% conducting component
sigma = sigma_w*(phi)^s; % S/m

% Compute skin depth from Schroeder et al 2015
del = sqrt(1/(pi*fc*sigma*mu_o)); %m

deltaAtten;
% Correct for portion of the signal likely due to scattering
deltaAtten_water = deltaAtten * (1-scattFrac);

deltaAtten_water = 10^(deltaAtten_water/10);

if deltaAtten_water == 0
    T_water = 0;
else
    T_water_iceMix = (del/2)*log(1/deltaAtten_water);
    T_water = T_water_iceMix*phi;
end




end

