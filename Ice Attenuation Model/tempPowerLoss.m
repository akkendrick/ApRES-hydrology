% Compute attenuation in dB from changing surface temperature
% Alex Kendrick
% 2/6/2017

% Compute attenuation rates as a function of temperature
[T_cel NgreenTot] = iceAttenuationModel();

 Tice = 50; %Ice thickness in m
 z = 0:0.1:Tice;

% Specify month for ice temperature model calculation
%months = [1,3,5,7,9,11];
%months = [1,2,3,4,5,6,7,8,9,10,11,12];

% Compute cumulative attenuation through top 50m over a timescale specified
% by fraction of months
months = 0:0.01:12;
for j = 1:length(months)
    month = months(j);
    [z tempDepthProfile tempTimeProfile] = iceTempModel(month, z);
    
    % Now combine both models to compute the attenuation rate as a func of
    % depth

    for i = 1:length(tempDepthProfile)
        attenProfile(i) = mean(NgreenTot(abs(T_cel - tempDepthProfile(i)) < 0.002));
    end
    
    allTempProfiles(j,:) = tempDepthProfile;
    allAttenProfiles(j,:) = attenProfile;
 
    delZ = z(2)-z(1); % Z spacing of temp profile in meters
    delZ = delZ * 10^-3; % Z spacing of temp profile in km
    cumulativeAtten = attenProfile*delZ;

%     figure(21)
%     plot(cumulativeAtten, z)
%     set(gca,'YDir','reverse')
%     
    % Two-way travel attenuation
    totalAtten(j) = 2*sum(cumulativeAtten);
    
end

%totalAtten

%%
close all
nSteps = size(allAttenProfiles,1);
colors = parula(13);
%months = [1,3,5,7,9,11];

% Plot the attenuation and temperature profiles through the top 50m of ice
% for each time step
% NOTE: this can take a long time depending on the number of timesteps
for i = 1:nSteps

    colorInd = floor(months(i))+1;
%     colorInd
%     colors(colorInd,:)
    
    figure(4)
    hold on
    plot(allAttenProfiles(i,:),z,'Color',colors(colorInd,:));
    set(gca,'YDir','reverse')
    xlabel('Attenuation rate N (dB/km, one way)')
    ylabel('Depth (m)')
    title('Attenuation Rate Profile')
    set(gca,'FontSize',14)
    legend('Jan', 'Mar','May','July','Sept','Nov')
    %legend('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

    figure(6)
    hold on
    plot(allTempProfiles(i,:),z,'Color',colors(colorInd,:));
    set(gca,'YDir','reverse')
    title('Temperature Depth Profile')
    xlabel('Temperature \circC')
    ylabel('Depth (m)')
    set(gca,'FontSize',14)
    legend('Jan', 'Mar','May','July','Sept','Nov')
    %legend('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

end
% colors = num2cell(winter(nSteps),2);
% set(h, {'color'},colors);


% colors = num2cell(winter(nSteps),2);
% set(h, {'color'},colors);


%%
% Plot the cumulative attenuation as a function of time
nSteps = size(allAttenProfiles,1);

% Set specifically for: %months = 0:0.01:12;
time3 = datetime(2014,1,1,1,1:438:525600,0);

figure(22)
hold on
plot(time3,totalAtten(1:end-1))
%plot(time3,totalAtten(1:end-1))
%plot(months,totalAtten)

xlabel('Month')
ylabel('Cumulative attenuation through top layer (dB)')
set(gca,'FontSize',14)
legend('Real Params','Modified Params')


% Compare to annual temp at weather station
load('AWS_tempData.mat')
date = AWS_date;
%date = cell2mat(date);
temp1 = AWS_temp1;
temp2 = AWS_temp2;

t = datetime(date,'InputFormat','yyyy-MM-dd');
%time2 = datetime(2014,1,1,1,1:438:525600,0);
% time3 = datetime(2014,1,1,1,1:438:525600,0);
% 
% years = [2014,2014,2014,2014,2014,2014];
% days = [1,1,1,1,1,1];
% % 
% % years = [2014,2014,2014,2014,2014,2014,2014,2014,2014,2014,2014,2014];
% % days = [1,1,1,1,1,1,1,1,1,1,1,1];
% 
% time2 = datetime(years,months,days);

figure(5)
hold on
plot(t,temp1,'.','MarkerSize',10)
plot(t,temp2,'.','MarkerSize',10)

plot(time3, tempTimeProfile(1:end-1),'m','LineWidth',2)
title('Annual Temperature at z = 0')
xlabel('Time (year fraction)')
ylabel('Temperature \circC')
set(gca,'FontSize',14)



% Quickly compute the temperature/atten difference in the top 15m from
% the mean value at 50m

tempDiff  = zeros(nSteps,151);
attenDiff = zeros(nSteps,151);


% Plot ice column temperature at a depth of 2m
twoMeterAtten = allAttenProfiles(:,z == 2);
twoMeterTemp = allTempProfiles(:,z == 2);


figure()
subplot(2,1,1)
plot(time3,twoMeterAtten(1:end-1))
xlabel('z (m)')
ylabel('Attenuation rate at 2m (dB/km)')

subplot(2,1,2)
plot(time3,twoMeterTemp(1:end-1))
xlabel('z (m)')
ylabel('Temperature at 2m (deg C)')

% Compute average ice temperature change over all depths?

tempDiff = zeros(1,length(z));
attenDiff = zeros(1,length(z));
for a = 1:length(z)
    currentTemp = allTempProfiles(:,z == z(a));
    currentAtten = allAttenProfiles(:,z==z(a));
    
    tempDiff(a) = max(currentTemp) - min(currentTemp);
    attenDiff(a) = max(currentAtten) - min(currentAtten);
    
end

disp('Average temperature difference in top 50m')
mean(tempDiff)

disp('Average attenuation difference in top 50m')
mean(attenDiff)


figure()
subplot(2,1,1)
plot(z,tempDiff)
xlabel('z (m)')
ylabel('Temperature Difference (deg C)')

subplot(2,1,2)
plot(z,attenDiff)
xlabel('z (m)')
ylabel('Attenuation difference (dB/km)')

% for i = 1:nSteps
%     
%     maxTemp = max(allTempProfiles(i,:));
%     meanTemp = allTempProfiles(i,500);
%     tempDiff(i) = maxTemp - meanTemp;
%     
%     maxAtten = max(allAttenProfiles(i,:));
%     meanAtten = allAttenProfiles(i,500);
%     attenDiff(i) = maxAtten - meanAtten;
% end
% 
% for i = 1:nSteps
%     for j = 1:151
%         meanTemp = allTempProfiles(i,500);
%         tempDiff(i,j) = allTempProfiles(i,j) - meanTemp;
%     
%         meanAtten = allAttenProfiles(i,500);
%         attenDiff(i,j) = allAttenProfiles(i,j) - meanAtten;
%     end
% end
% 
% meanTempDiff = mean(tempDiff,2);
% meanAttenDiff = mean(attenDiff,2);

