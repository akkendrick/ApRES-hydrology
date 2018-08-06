function scattFrac = computeScattFrac(phi, runoff, sigma_w, r)
% Compute fraction of attenuation attributable to volumetric scattering
%  
%	Inputs
%		phi - porosity of water storage
%		runoff - runoff fraction
%		sigma_w - electrical conductivity of water filled pores
%		r - pore radius
%
%	Outputs
%		scattFrac - fraction of attenuation attributable to volumetric scattering


f = 300*10^6; %Hz

startPower = -40.00; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modeledMelt = 0.01:0.01:10;

for a = 1:length(modeledMelt)
    [waterAtten, sigma_eff, T] = waterIceMixReturn(modeledMelt(a),f,sigma_w,phi,runoff);
    [Qsca, scatteringAtten] = computeMie(r,T,phi); 
    
    modeledAtten(a+1) = startPower + scatteringAtten + waterAtten;
    
    scattFrac(a+1) = scatteringAtten / (scatteringAtten+waterAtten);
end
% 
% figure(1)
% hold on
% plot(modeledMelt, scattFrac(2:end))

scattFrac = scattFrac(2);
end

