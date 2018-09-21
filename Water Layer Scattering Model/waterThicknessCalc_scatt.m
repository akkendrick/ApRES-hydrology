% Estimate water layer thickness for a given attenuation in dB

% Define ApRES and porous layer parameters
fc = 300*10^6;
sigma_w = 0.00016;
phi = 0.3;
r = 0.034;
runoff = 0;

% Compute fraction of attenuation attributable to volumetric scattering
scattFrac = computeScattFrac(phi,runoff,sigma_w,r);

% Specify change in attenuation due to proposed water storage
deltaAtten = (80-28)*-1.0;

% Compute water + porous layer thickness from specifed change in attenuation
[T_water, T_water_iceMix] = estimateMeltScatt(deltaAtten,scattFrac,fc,sigma_w,phi);

disp('For an attenuation of:')
deltaAtten
disp('We estimate water thickness of:')
T_water
disp('And a water storage region thickness of:')
T_water_iceMix