function [ z tempDepthProfile tempTimeProfile ] = iceTempModel(month, z)
% Compute surface ice temperature based on derivation in Hooke 2005
% Alex Kendrick
% 1/30/2017

% Inputs
%   month - number corresponding to the month of a given year (ie 1 for
%   Jan)
% Outputs
%   tempDepthProfile - temperature as a function of the first 15m of ice


% Implementing surface temperature model from Hooke, Principles of Glacier
% Mechanics 2005


kappa = 37.2; % m^2/a

% Modify kappa drastically for unrealistic worst case scenario
%kappa = 300;

meanTemp = -6.54;

%annualTempAmp = 16.79; %computed by looking at half an annual cycle of Greenland temp data

% Modify annualTempAmp for worst case scenario
%annualTempAmp = 20;
annualTempAmp = 40;

omega = 2.8*pi; %2 pi radians per year


z0 = 0;
timeOrig = 0:0.01:12;

boxFunc = ones(1, length(timeOrig));
boxFunc(timeOrig > 11.1) = 0;
boxFunc(timeOrig < 2.6) = 0;

time = timeOrig/12;


tempTimeProfile = meanTemp + 1/2.*(annualTempAmp.*exp(-z0.*sqrt(omega/(2*kappa))).*sin(omega.*(time-0.25)...
    - z0.*sqrt(omega/(2.*kappa))-1.));
tempTimeProfile = tempTimeProfile .* boxFunc;

%tempTimeProfile(timeOrig > 10.5) = -14.74;
%tempTimeProfile(timeOrig > 10.5) = -16.31;
tempTimeProfile(timeOrig > 10.5) = -26.08;



%tempTimeProfile(timeOrig < 2.6) = -14.61;
%tempTimeProfile(timeOrig < 2.6) = -16.16;
tempTimeProfile(timeOrig < 2.6) = -25.77;


t = month;
if t > 10.5
    t = 10.5;
end

if t < 2.6
    t = 2.6;
end

t = t/12;


tempDepthProfile = meanTemp+1/2.*(annualTempAmp.*exp(-z.*sqrt(omega/(2*kappa))).*sin(omega*(t-0.25)...
- z.*sqrt(omega/(2.*kappa))));

% figure(4)
% subplot(1,2,1)
% hold on
% plot(tempDepthProfile,z)
% set(gca,'YDir','reverse')
%     
% figure(4)
% subplot(1,2,1)
% title('Temperature Depth Profile')
% xlabel('Temperature \circC')
% ylabel('Depth (m)')
% legend('Jan','Mar','May','July','Sept','Nov')
% set(gca,'FontSize',14)
% 
% 
% figure(5)
% hold on
% plot(time, tempTimeProfile)
% title('Annual Temperature at z = 0')
% xlabel('Time (year fraction)')
% ylabel('Temperature \circC')
% set(gca,'FontSize',14)


end