% plotRadarAttenuation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to plot the observed ApRES radar attenuation, adjusting for a
% constant offset between deployment 1 and deployment 2
%
% Alex Kendrick
% 7/27/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('radarData_deploy1_reprocess.mat')
meanRadarDat_deploy1 = meanRadarDat;
t_deploy1 = t(1:end,1);
yReflector_deploy1 = yReflector;
xLoc_deploy1 = xLoc;
yLoc_deploy1 = yLoc;

load('radarData_deploy2_reprocess.mat')
meanRadarDat_deploy2 = meanRadarDat;
t_deploy2 = t(1:end,1);
yReflector_deploy2 = yReflector;
xLoc_deploy2 = xLoc;
yLoc_deploy2 = yLoc;

% Come up with linear scaling relation between deployment 1 and deployment
% 2
% Identify linear portion at end of deployment 1 and start of deployment 2
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
newDates = {'2014-07-22','2014-07-24','2014-07-26','2014-07-28','2014-07-30',...
    '2014-08-02'};
newTimeVec = datetime(newDates,'InputFormat','yyyy-MM-dd');

% Compute estimated data for entire time window with average fitting
% parameters
combinedTimeVec = [tSubset1; newTimeVec'; tSubset2];
combinedTimeVecNum = datenum(combinedTimeVec);
combinedFit = averageFit(1)*combinedTimeVecNum+averageFit(2);

% Adjust starting position in line with deployment 1 data
scaleFactor = fittedY(1) - combinedFit(1);
combinedFit = combinedFit + scaleFactor;

% Compute average offset needed to shift deployment 2
rawOffset = combinedFit(309:end) - meanRadarDat_deploy2(1,1:98)';
offsetCorr = mean(rawOffset)

figure(8)
grid on
box
title('Offset correction between deployment 1 and 2')
hold on
plot(t_deploy2, meanRadarDat_deploy2,'.')

% Now scale deployment2 to match modelled data
% For reprocessed data correction is: 38.4925 dB
meanRadarDat_deploy2 = meanRadarDat_deploy2 + offsetCorr;
yReflector_deploy2 = yReflector + offsetCorr;

figure(8)
plot(t_deploy1, meanRadarDat_deploy1,'.')
plot(t_deploy2, meanRadarDat_deploy2,'.')
plot(combinedTimeVec, combinedFit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load bed power and apply same correction 
load('radarData_deploy1_bed_reprocess.mat')
yLoc_bed = 617;
xLoc_bed = -2.87;

bedPower_deploy1 = baseYReflector;

load('radarData_deploy2_bed_reprocess.mat')
bedPower_deploy2 = baseYReflector;

% Now scale deploy2 bed power to match 
bedPower_deploy2 = bedPower_deploy2 + offsetCorr;

figure(20)
grid on
box
title('Bed power, average internal layer comparision')
hold on
plot(t_deploy1, bedPower_deploy1,'.r','MarkerSize',10)
plot(t_deploy1, meanRadarDat_deploy1,'.b','MarkerSize',10)

plot(t_deploy2, bedPower_deploy2,'.r','MarkerSize',10)
plot(t_deploy2, meanRadarDat_deploy2,'.b','MarkerSize',10)
ylim([-120, -20])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot internal layers
figure(9)
grid on
box
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

% create the colormap:
numcolors = length(yLoc_deploy2);
% create the gradients
clrmap = cell2mat(arrayfun(@(a,b)linspace(a,b,numcolors )',color1,color2,'uni',false));

for j = 1:length(yLoc_deploy2)
       
    figure(9)
    hold on
    title('dB vs Time of Reflectors in the Y-Profile')
    xlabel('Time')
    ylabel('dB (Vms)')
    plot(t_deploy2, yReflector_deploy2(j,:), '.','Color',clrmap(j,:))
    set(gca, 'FontSize',16)
    ylim([-95,-35])
    
end

plot(t_deploy1, meanRadarDat_deploy1,'.b','MarkerSize',10)
plot(t_deploy2, meanRadarDat_deploy2,'.b','MarkerSize',10)

% Plot fitting for reference
% plot(t_deploy1(1335:end),fittedY)
% plot(t_deploy2(1:98),fittedY2)
% plot(combinedTimeVec, combinedFit,'k','LineWidth',3)

alpha(0.5)
ylim([-110,-35])
xlim([datetime(2014,05,01),datetime(2014,12,01)])
set(gca,'FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


