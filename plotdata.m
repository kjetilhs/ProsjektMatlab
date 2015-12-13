close all;
clear;

PostPro = load('log_folder/PostPro/rtklib_ubxstream_x8_log201512081128.pos');
flightdata = load('log_folder/rtklibLog/rtklib_output201512081128.pos');
len_msg = length(flightdata(:,1));

tStart = 1542;
tStop = 4958;

t = flightdata(tStart:tStop,2)-flightdata(tStart,2);
e = flightdata(tStart:tStop,3);
n = flightdata(tStart:tStop,4);
d = -flightdata(tStart:tStop,5);
quality = flightdata(tStart:tStop,6);
num_sat = flightdata(tStart:tStop,7);
sdx = flightdata(tStart:tStop,8);
sdy = flightdata(tStart:tStop,9);
sdz = flightdata(tStart:tStop,10);
sdxy = flightdata(tStart:tStop,11);
sdyz = flightdata(tStart:tStop,12);
sdzx = flightdata(tStart:tStop,13);
age = flightdata(tStart:tStop,14);
ratio = flightdata(tStart:tStop,15);
ve = flightdata(tStart:tStop,16);
vn = flightdata(tStart:tStop,17);
vd = -flightdata(tStart:tStop,18);

tp = PostPro(tStart:tStop,2)-PostPro(tStart,2);
ep = PostPro(tStart:tStop,3);
np = PostPro(tStart:tStop,4);
dp = -PostPro(tStart:tStop,5);
qualityp = PostPro(tStart:tStop,6);
num_satp = PostPro(tStart:tStop,7);
sdxp = PostPro(tStart:tStop,8);
sdyp = PostPro(tStart:tStop,9);
sdzp = PostPro(tStart:tStop,10);
sdxyp = PostPro(tStart:tStop,11);
sdyzp = PostPro(tStart:tStop,12);
sdzxp = PostPro(tStart:tStop,13);
agep = PostPro(tStart:tStop,14);
ratiop = PostPro(tStart:tStop,15);
%% Fixing timestamp in tow
start = 1;
stop = 0;
for i=1:length(t)-1
    if t(i) == t(start)
        stop = stop + 1;
    elseif t(i) ~= t(start)
        
        diff = i-start;
        if diff == 0
            diff = 1;
        end
        for k=0:diff-1
           t(k+start) = t(k+start)+k/diff; 
        end
        start = i;
        stop = i;
    end
end
errord = d-dp;
errore = e-dp;
errorn = n-np;
%% Retrive FIX,FLOAT and NONE solution
numFix = flightdata(tStart:tStop,6) == 1;
[rowStart,~] = find(numFix);
numFix = length(rowStart);
tFix = zeros(1,numFix);
eFix = zeros(1,numFix);
nFix = zeros(1,numFix);
dFix = zeros(1,numFix);
qualityFix = zeros(1,numFix);
num_satFix = zeros(1,numFix);
sdxFix = zeros(1,numFix);
sdyFix = zeros(1,numFix);
sdzFix = zeros(1,numFix);
sdxyFix = zeros(1,numFix);
sdyzFix = zeros(1,numFix);
sdzxFix = zeros(1,numFix);
ageFix = zeros(1,numFix);
ratioFix = zeros(1,numFix);
veFix = zeros(1,numFix);
vnFix = zeros(1,numFix);
vdFix = zeros(1,numFix);


numFloat = flightdata(tStart:tStop,6) == 1;
[rowStart,~] = find(numFloat);
numFloat = length(rowStart);
tFloat = zeros(1,numFloat);
eFloat = zeros(1,numFloat);
nFloat = zeros(1,numFloat);
dFloat = zeros(1,numFloat);
qualityFloat = zeros(1,numFloat);
num_satFloat = zeros(1,numFloat);
sdxFloat = zeros(1,numFloat);
sdyFloat = zeros(1,numFloat);
sdzFloat = zeros(1,numFloat);
sdxyFloat = zeros(1,numFloat);
sdyzFloat = zeros(1,numFloat);
sdzxFloat = zeros(1,numFloat);
ageFloat = zeros(1,numFloat);
ratioFloat = zeros(1,numFloat);
veFloat = zeros(1,numFloat);
vnFloat = zeros(1,numFloat);
vdFloat = zeros(1,numFloat);

numNone = tStop-tStart-numFix-numFloat;
if numNone <0
    numNone = 0;
end

tNone = zeros(1,numNone);
eNone = zeros(1,numNone);
nNone = zeros(1,numNone);
dNone = zeros(1,numNone);
qualityNone = zeros(1,numNone);
num_satNone = zeros(1,numNone);
sdxNone = zeros(1,numNone);
sdyNone = zeros(1,numNone);
sdzNone = zeros(1,numNone);
sdxyNone = zeros(1,numNone);
sdyzNone = zeros(1,numNone);
sdzxNone = zeros(1,numNone);
ageNone = zeros(1,numNone);
ratioNone = zeros(1,numNone);
veNone = zeros(1,numNone);
vnNone = zeros(1,numNone);
vdNone = zeros(1,numNone);
fix = 1;
float = 1;
none = 1;
for i=1:tStop-tStart
   if quality(i) == 1 %FIX
       quality(i) = 3;
       qualityFix(fix) = 3;
       tFix(fix) = t(i);
       eFix(fix) = e(i);
       nFix(fix) = n(i);
       dFix(fix) = d(i);
       num_satFix(fix) = num_sat(i);
       sdxFix(fix) = sdx(i);
       sdyFix(fix) = sdy(i);
       sdzFix(fix) = sdz(i);
       sdxyFix(fix) = sdxy(i);
       sdyzFix(fix) = sdyz(i);
       sdzxFix(fix) = sdzx(i);
       ageFix(fix) = age(i);
       ratioFix(fix) = ratio(i);
       veFix(fix) = ve(i);
       vnFix(fix) = vn(i);
       vdFix(fix) = vd(i);
       fix = fix + 1;
   elseif quality(i) == 2 %FLOAT
       quality(i) = 2;
       qualityFloat(float) = 2;
       tFloat(float) = t(i);
       eFloat(float) = e(i);
       nFloat(float) = n(i);
       dFloat(float) = d(i);
       num_satFloat(float) = num_sat(i);
       sdxFloat(float) = sdx(i);
       sdyFloat(float) = sdy(i);
       sdzFloat(float) = sdz(i);
       sdxyFloat(float) = sdxy(i);
       sdyzFloat(float) = sdyz(i);
       sdzxFloat(float) = sdzx(i);
       ageFloat(float) = age(i);
       ratioFloat(float) = ratio(i);
       veFloat(float) = ve(i);
       vnFloat(float) = vn(i);
       vdFloat(float) = vd(i);
       float = float + 1;
   else %NONE
       quality(i) = 0;
       qualityNone(none) = 0;
       tNone(none) = t(i);
       eNone(none) = e(i);
       nNone(none) = n(i);
       dNone(none) = d(i);
       num_satNone(none) = num_sat(i);
       sdxNone(none) = sdx(i);
       sdyNone(none) = sdy(i);
       sdzNone(none) = sdz(i);
       sdxyNone(none) = sdxy(i);
       sdyzNone(none) = sdyz(i);
       sdzxNone(none) = sdzx(i);
       ageNone(none) = age(i);
       ratioNone(none) = ratio(i);
       veNone(none) = ve(i);
       vnNone(none) = vn(i);
       vdNone(none) = vd(i);
       none = none + 1;
   end
   if qualityp(i) ==1 %Fix
       qualityp(i) = 3;
   elseif qualityp(i) == 2 %Float
       qualityp(i) = 2;
   else %None
       qualityp(i) = 0;
   end
end
xyStart = 2500;
xyStop = 3000;
%% Plots
figure(1);
    plot3(e(xyStart:xyStop),n(xyStart:xyStop),d(xyStart:xyStop),'xb');
%     hold on;
%     plot(eFloat(xyStart:xyStop),nFloat(xyStart:xyStop),'xr');
%     plot(eNone,nNone,'xg');
    grid on;
    title('XY'); 
    xlabel('East [m]');
    ylabel('North [m]');
    zlabel('Down [m]');
%     legend('Fix','Float','None');

    
    figure(2);
    subplot(2,1,1);
    plot(tFix,dFix,'xb');
    hold on;
    plot(tFloat,dFloat,'xr');
    plot(tNone,dNone,'xg');
    plot(tp,dp,'--c');
    grid on;
    title('Down');
    legend('Fix','Float','None','Post processed')
    xlabel('Time [s]');
    ylabel('Meter [m]');

    subplot(2,1,2)
    plot(t,quality,'b');
    hold on;
    plot(tp,qualityp,'r');
%     plot(PostTime-PostTime(1),PostType,'--g');
    grid on;
    title('Integer ambiguity solution: Fix = 3, Float = 2, None = 0')
    ylabel('Solution type')
    xlabel('Time [s]');
    legend('Real time','Post processed');
    ylim([0 5]);

    figure(3)
    subplot(3,1,1)
    plot(tFix,vnFix,'xb');
    hold on;
    plot(tFloat,vnFloat,'xr');
    plot(tNone,vnNone,'xg');
    grid on;
    legend('Fix','Float','None');
    title('Velocity in North direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,2)
    plot(tFix,veFix,'xb');
    hold on;
    plot(tFloat,veFloat,'xr');
    plot(tNone,veNone,'xg');
    grid on;
    legend('Fix','Float','None');
    title('Velocity in East direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,3)
    plot(tFix,vdFix,'xb');
    hold on;
    plot(tFloat,vdFloat,'xr');
    plot(tNone,vdNone,'xg');
    grid on;
    legend('Fix','Float','None');
    title('Velocity in Down direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    figure(4)
    plot(t,num_sat,'b');
    hold on;
%     plot(timestamp_r(1:TimeEndr)-timeStart,satellites_r(1:TimeEndr),'r');
    grid on;
    title('Visable satellites')
    xlabel('Time [s]');
    ylabel('Number of satellites');
%     legend('Rtklib');
    ylim([4 12])
%     
%     %% Histogram
%     figure(5)
%     subplot(2,1,1);
%     histogram(deltatimetow_r)
%     title('RTKLIB tow');
%     xlabel('Time [s]');
%     
%     subplot(2,1,2);
%     histogram(deltatime_r)
%     title('RTKLIB timestamp');
%     xlabel('Time [s]');

figure(6)
subplot(2,1,1);
plot(tp,errord);
grid on;
title('Error in Down position');
xlabel('Time [s]');
ylabel('Meter [m]');
ylim([-2 2]);
subplot(2,1,2)
plot(t,quality,'b');
hold on;
plot(tp,qualityp,'r');
grid on;
title('Integer ambiguity solution: Fix = 3, Float = 2,Obs = 1, None = 0')
ylabel('Solution type')
xlabel('Time [s]');
legend('Real time','Post processed');
ylim([0 5]);

figure(7)
subplot(3,1,1);
plot(tp,errord);
grid on;
title('Error in North position');
xlabel('Time [s]');
ylabel('Meter [m]');
ylim([-2 2]);
subplot(3,1,2);
plot(tp,errord);
grid on;
title('Error in East position');
xlabel('Time [s]');
ylabel('Meter [m]');
ylim([-2 2]);
subplot(3,1,3);
plot(tp,errord);
grid on;
title('Error in Down position');
xlabel('Time [s]');
ylabel('Meter [m]');
ylim([-2 2]);

figure(8);
subplot(2,1,1);
plot(tFix,nFix,'xb');
hold on;
plot(tFloat,nFloat,'xr');
plot(tNone,nNone,'xg');
plot(tp,np,'--c');
grid on;
title('North');
legend('Fix','Float','None','Post processed')
xlabel('Time [s]');
ylabel('Meter [m]');
subplot(2,1,2);
plot(tFix,eFix,'xb');
hold on;
plot(tFloat,eFloat,'xr');
plot(tNone,eNone,'xg');
plot(tp,ep,'--c');
grid on;
title('East');
legend('Fix','Float','None','Post processed')
xlabel('Time [s]');
ylabel('Meter [m]');