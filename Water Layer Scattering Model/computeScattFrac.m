function scattFrac = computeScattFrac(phi, runoff, sigma_w, r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

eps_i = 3.17;
eps_w = 80;
eps_b = 18;
f = 300*10^6; %Hz

startPower = -40.00; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modeledMelt = 0.01:0.01:10;

for a = 1:length(modeledMelt)
    [waterAtten, sigma_eff, T] = waterIceMixReturn(modeledMelt(a),f,sigma_w,eps_w,eps_b,eps_i,phi,runoff);
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

