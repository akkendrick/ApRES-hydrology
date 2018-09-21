% modelEnglacialAttenuation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to plot the observed ApRES radar attenuation, and model the
% estimated attenuation expected from the 2014 AWS melt data. The script
% estimates attenuation for a range of pore sizes given a porosity 
% and conductivity of the englacial water. The best fit forward model is
% determined.
%
% Alex Kendrick
% 7/27/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift by original offset to measure relative power changes in dB
initialOffset = 40;

load('AWS_meltData.mat')

load('internalLayer_radarData_deploy1.mat')
meanRadarDat_deploy1 = meanRadarDat + initialOffset;
t_deploy1 = t(1:end,1);
yReflector_deploy1 = yReflector + initialOffset;
xLoc_deploy1 = xLoc;
yLoc_deploy1 = yLoc;

load('internalLayer_radarData_deploy2.mat')
meanRadarDat_deploy2 = meanRadarDat + initialOffset;
t_deploy2 = t(1:end,1);
yReflector_deploy2 = yReflector + initialOffset;
xLoc_deploy2 = xLoc;
yLoc_deploy2 = yLoc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify porosity for water/ice mixture
% Now make an array of phi for a given conductivity to perform best fit
phi = 0.175;
runoff = 0;
T_sum_AWS = 0;

%sigma_w = 0.007;
sigma_w = 0.0003;

% Loop over a given set of pore radii 
r  = 0.034:0.001:0.043;
f = 300*10^6; %Hz

startPower = 0.00; 

eps_i = 3.17;
eps_w = 80;
eps_b = 18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


newTimeVecNum = datenum(newTimeVec);
projectedData = averageFit(1)*newTimeVecNum + averageFit(2);
projectedData = projectedData';

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
yReflector_deploy2 = yReflector_deploy2 + offsetCorr;

figure(8)
plot(t_deploy1, meanRadarDat_deploy1,'.')
plot(t_deploy2, meanRadarDat_deploy2,'.')
plot(combinedTimeVec, combinedFit)

totalTimeVec = [t_deploy1; newTimeVec'; t_deploy2];
totalMeanRadarAtten = [meanRadarDat_deploy1 projectedData' meanRadarDat_deploy2];

% Interp attenuation data onto AWS dates
totalAWSTimeVecNum = datenum(totalTimeVec);
AWSdateNums = datenum(dailyAWSDate);
meanAWSAttenInterp = interp1(totalAWSTimeVecNum,totalMeanRadarAtten,AWSdateNums);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load bed power and apply same correction 
load('bed_radarData_deploy1.mat')
yLoc_bed = 617;
xLoc_bed = -2.87;

bedPower_deploy1 = baseYReflector + initialOffset;

load('bed_radarData_deploy2.mat')
bedPower_deploy2 = baseYReflector + initialOffset;

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
ylim([-70, 15])


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
    
end

plot(t_deploy1, meanRadarDat_deploy1,'.b','MarkerSize',10)
plot(t_deploy2, meanRadarDat_deploy2,'.b','MarkerSize',10)

% Plot fitting for reference
% plot(t_deploy1(1335:end),fittedY)
% plot(t_deploy2(1:98),fittedY2)
% plot(combinedTimeVec, combinedFit,'k','LineWidth',3)

alpha(0.5)
ylim([-70, 10])
xlim([datetime(2014,05,01),datetime(2014,12,01)])
set(gca,'FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate attenuation from storage
meltAWS = zeros(length(dailyAWSDate)+1,length(r));
totalAblationAWS = zeros(length(dailyAWSDate),length(r));
attenModelAWS = zeros(length(dailyAWSDate),length(r));
AWS_ablation = zeros(length(dailyAWSDate),length(r));

% Loop over timeVec and calculate ablation for each pore size
for b = 1:length(r)
    attenModelAWS(1,b) = startPower;
    T_sum_RACMO = 0;

    for j=1:length(dailyAWSDate)


        meltAWS(j) = dailyAWSMelt(j); %
        totalAblationAWS(j+1,b) = totalAblationAWS(j,b) + meltAWS(j);
        
        [waterAtten, sigma_eff, T] = waterIceMixReturn(totalAblationAWS(j),f,sigma_w,phi,runoff);
        [Qsa, scatteringAtten] = computeMie(r(b),T,phi);
        
        attenModelAWS(j+1,b) = startPower + waterAtten + scatteringAtten;

        AWS_ablation(j+1,b) = T ;
       if ~isnan(T)
           T_sum_AWS = T;
       end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute error between AWS attenuation curve and observed mean internal
% layer data

goodAWSInterp = ~isnan(meanAWSAttenInterp);
goodAWSInterpDates = dailyAWSDate(goodAWSInterp);
meanAWSAttenInterp = meanAWSAttenInterp(goodAWSInterp);
goodAttenModelAWS = attenModelAWS(goodAWSInterp,:);

AWSRMSError = zeros(length(goodAWSInterpDates),length(phi));

for b = 1:length(r)
    AWSRMSError(:,b) = sqrt((goodAttenModelAWS(:,b) - meanAWSAttenInterp).^2);
end

summedAWSError = sum(AWSRMSError);

[minAWSError, AWSindex] = min(summedAWSError);
bestAWSFitCurve = goodAttenModelAWS(:, AWSindex);

fprintf('AWS Best fit model for sigma %5.4f\n',sigma_w)
fprintf('AWS Porosity of %5.3f\n',phi)
fprintf('AWS Pore Size of %5.3f\n',r(AWSindex))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxColor = [100 170 45]/255;
bestColor = [37 37 37]/255;
minColor = [40 235 200]/255;
meanColor = [165 15 21]/255;
RACMOColor = [205 145 60]/255;


figure(9)
% Plot AWS Results
h1d = plot(goodAWSInterpDates, goodAttenModelAWS(:,1),'Color',minColor,'LineWidth',5); h1d.Color(4) = 0.75; % min porosity
h1c = plot(goodAWSInterpDates, goodAttenModelAWS(:,end),'Color',maxColor,'LineWidth',5); h1c.Color(4) = 0.75; % max porosity
h1b = plot(goodAWSInterpDates, bestAWSFitCurve,'Color',bestColor,'LineWidth',5); h1b.Color(4) = 0.75;


