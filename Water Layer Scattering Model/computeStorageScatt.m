% Estimate water storage from observed internal layer power drop

% Compute estimated attenuation curve
startPower = -40.00; 
eps_w = 80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify porosity for water/ice mixture
phi = 0.3;
SMBrunoff = 0;
T_sum_AWS = 0;
sigma_w = 0.00016;
%sigma_w = 0.00399;
fc = 300*10^6; %Hz
r = 0.04;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load internal layers
load('radarData_deploy1_reprocess.mat')
meanRadarDat_deploy1 = meanRadarDat;
t_deploy1 = t;
yReflector_deploy1 = yReflector;
xLoc_deploy1 = xLoc;
yLoc_deploy1 = yLoc;

load('radarData_deploy2_reprocess.mat')
meanRadarDat_deploy2 = meanRadarDat;
t_deploy2 = t(:,1);
yReflector_deploy2 = yReflector;
xLoc_deploy2 = xLoc;
yLoc_deploy2 = yLoc;

% Now scale deployment2 to match modelled data
meanRadarDat_deploy2 = meanRadarDat_deploy2 - 38.4925;
yReflector_deploy2 = yReflector - 38.4925;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot internal layer power  
figure(9)
hold on

% Plot the internal layers coloring based on depth

color1 = [239 243 255]/255;
color2 = [8 81 156]/255;

numcolors = length(yLoc_deploy1);
% create the gradients
clrmap = cell2mat(arrayfun(@(a,b)linspace(a,b,numcolors )',color1,color2,'uni',false));

for j = 1:length(yLoc_deploy1)
       
    figure(9)
    hold on

    plot(t_deploy1, yReflector_deploy1(j,:), '.','Color',clrmap(j,:))
    set(gca, 'FontSize',14)
  
end

numcolors = length(yLoc_deploy2);
% create the gradients
clrmap = cell2mat(arrayfun(@(a,b)linspace(a,b,numcolors )',color1,color2,'uni',false));

for j = 1:length(yLoc_deploy2)
       
    figure(9)
    hold on
    box
    grid on
    title('dB vs Time of Reflectors in the Y-Profile')
    xlabel('Time')
    ylabel('dB (Vms)')
    plot(t_deploy2, yReflector_deploy2(j,:), '.','Color',clrmap(j,:))
    ylim([-95,-35])
  
end

plot(t_deploy1, meanRadarDat_deploy1,'.b','MarkerSize',10)
plot(t_deploy2, meanRadarDat_deploy2,'.b','MarkerSize',10)

alpha(0.5)
ylim([-110,-35])
xlim([datetime(2014,05,01),datetime(2014,12,01)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First stitch together the data sets
tSubset1 = t_deploy1(1335:end);
tSubsetNum1 = datenum(tSubset1);

tSubset2 = t_deploy2(1:98);
tSubsetNum2 = datenum(tSubset2);

ySubset = meanRadarDat_deploy1(1335:end);
ySubset2 = meanRadarDat_deploy2(1:98);

% Fit the data
p1 = polyfit(tSubsetNum1, ySubset', 1);
p2 = polyfit(tSubsetNum2, ySubset2',1);

% Trend for end of deployment 1
fittedY = p1(1)*tSubsetNum1+p1(2);

% Trend for start of deployment 2
fittedY2 = p2(1)*tSubsetNum2+p2(2);

% Average both trends
averageFit = (p1 + p2) ./ 2;

% Create time vector to account for gap in data
time3 = datetime(2014,7,17,1,1:1440:25920,0);
%time3 = datetime(2014,7,10,1,1:1440:43200,0);

newDates = time3;
newTimeVec = datetime(newDates,'InputFormat','yyyy-MM-dd');

newTimeVecNum = datenum(newTimeVec);
projectedData = averageFit(1)*newTimeVecNum + averageFit(2);

% Compute estimated data for entire time window with average fitting
% parameters
totalTimeVec = [t_deploy1; newTimeVec'; t_deploy2];
totalMeanRadarAtten = [meanRadarDat_deploy1 projectedData meanRadarDat_deploy2];

combinedTimeVecNum = datenum(totalTimeVec);
combinedFit = averageFit(1)*combinedTimeVecNum+averageFit(2);

% Adjust starting position in line with deployment 1 data
scaleFactor = fittedY(1) - combinedFit(1);
combinedFit = combinedFit + scaleFactor;

% Compute average offset needed to shift deployment 2
offsetCorr = 38.4925;

% Now scale deployment2 to match modelled data
% For reprocessed data correction is: 38.4925 dB
meanRadarDat_deploy2 = meanRadarDat_deploy2 + offsetCorr;
yReflector_deploy2 = yReflector + offsetCorr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now that we have the observed data, compute estimated storage from the
% combined observed attenuation

estimatedTotalMelt = zeros(length(totalTimeVec)+1,1);
predictedMelt = zeros(length(totalTimeVec)+1,1);
estimatedFirnThickness = zeros(length(totalTimeVec)+1,1);

scattFrac = computeScattFrac(phi, SMBrunoff, sigma_w, r);

powerDiffVec = zeros(length(totalTimeVec),1);
for a = 1:length(totalTimeVec)

    if a == 1
        prevPower = startPower;
    else
        prevPower = totalMeanRadarAtten(a-1);
    end
    
    powerDiff = startPower - totalMeanRadarAtten(a);
    powerDiffVec(a) = powerDiff;
    

    % Find the scattering fraction closest to the observed powerDiff
    [T_water, T_water_iceMix] = estimateMeltScatt(powerDiff,scattFrac,fc,sigma_w,phi);
    T_water = -1.0*T_water;
    T_water_iceMix = -1.0*T_water_iceMix;
    
    estimatedTotalMelt(a+1) = T_water;
    predictedMelt(a) = T_water;
    estimatedFirnThickness(a+1) = T_water_iceMix;
    
end


%Average runoff fraction
disp('For phi of')
phi

disp('Final water layer thickness')
estimatedTotalMelt(end-1)

disp('Final saturated firn layer thickness')
estimatedFirnThickness(end-1)

maxColor = [100 170 45]/255;
bestColor = [37 37 37]/255;
minColor = [40 235 200]/255;
meanColor = [165 15 21]/255;


figure(10)
subplot(2,1,1)
hold on
box
grid on
plot(totalTimeVec,(estimatedTotalMelt(1:end-1)),'.')
title('Estimated stored melt')

figure(10)
subplot(2,1,2)
hold on
box
grid on
plot(totalTimeVec, estimatedFirnThickness(1:end-1),'.');
title('Estimated firn thickness')
xlim([datetime('May 1, 2014'),datetime('Nov 30, 2014')])

save('predictedMelt.mat','totalTimeVec', 'estimatedTotalMelt','predictedMelt')
