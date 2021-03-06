 function [ iceAtten, sigma, T ] = waterIceMixReturn( T_water, fc, sigma_w,phi,runoff)
% Computes attenuation due to water/ice mix using a mixing model and
% attenuation through a conductor. See Geldsetzer et al 2009 for mixing
% model.
%
%	Inputs
%		T_water - water layer thickness
%		fc - center frequency
%		sigma_w - electrical conductivity of water
%       phi - porosity of water/ice mix
%       runoff - surface runoff fraction
%
%  Outputs
%       iceAtten - modeled attenuation through water/ice mixture
%       sigma - conductivity of water/ice mixture via mixing model
%       T - thickness of water/ice mixture


mu_o = 1.2566370614*10^(-6); %H/m;

if runoff > 0
    T_water = T_water * (1-runoff);
end

% Implement mixing model
% So if water is of thickness T, we need
T_water_iceMix = T_water/phi;%

% Skin depth matches results in Dowdeswell and Evans 2010
% Conductivity mixing model and mean value of s from Geldsetzer et al 2009
s = 1.67;

% The water ice mixture (phi or v) is considered to be the thickness of the
% conducting component
sigma = sigma_w*(phi)^s; % S/m

% Compute skin depth from Schroeder et al 2015
%disp(T)
del = sqrt(1/(pi*fc*sigma*mu_o)); %m

% Rb = abs((sqrt(eps_b)-sqrt(eps_w))/(sqrt(eps_b)+sqrt(eps_w)));
% Rw = abs((sqrt(eps_w)-sqrt(eps_i))/(sqrt(eps_w)+sqrt(eps_i)));

%waterAtten = Rw*exp(-2*T/del);
iceAtten = exp(-2*T_water_iceMix/del);
iceAtten = 10*log10(iceAtten);

% Force water layer attenuation to zero if T = 0
if T_water == 0
    iceAtten = 0;
end

% Return the thickness of the water/ice mixture
T = T_water_iceMix;

end
