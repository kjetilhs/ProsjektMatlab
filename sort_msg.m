%piksi 22
%rtklib 23
close all;
clear;
load('log_folder/neptusLog/Vandring6/Data.mat')
PostPro = load('log_folder/PostPro/rtklib_ubxstream_x8_log201512061145.pos');
PIXI = 1;
len = length(RtkFix.src_ent);

numberRTK = 0;
numberPIXI = 0;
%% Finding the number of RKT and PIXI
for i=1:len
    if RtkFix.src_ent(i) == 22
        numberPIXI = numberPIXI+1;
    elseif RtkFix.src_ent(i) == 23
        numberRTK = numberRTK +1;
    end
end
%% Initializing variables for RTK and PIXI

% RTK
timestamp_r = zeros(1,numberRTK);
src_r = zeros(1,numberRTK);
src_ent_r = zeros(1,numberRTK);
dst_r =zeros(1,numberRTK);
tow_r = zeros(1,numberRTK);
n_r = zeros(1,numberRTK);
e_r = zeros(1,numberRTK);
d_r = zeros(1,numberRTK);
v_n_r = zeros(1,numberRTK);
v_e_r = zeros(1,numberRTK);
v_d_r = zeros(1,numberRTK);
satellites_r = zeros(1,numberRTK);
iar_hyp_r = zeros(1,numberRTK);
iar_ratio_r = zeros(1,numberRTK);
type_r = zeros(1,numberRTK);

% PIXI
timestamp_p = zeros(1,numberPIXI);
src_p = zeros(1,numberPIXI);
src_ent_p = zeros(1,numberPIXI);
dst_p =zeros(1,numberPIXI);
tow_p = zeros(1,numberPIXI);
n_p = zeros(1,numberPIXI);
e_p = zeros(1,numberPIXI);
d_p = zeros(1,numberPIXI);
v_n_p = zeros(1,numberPIXI);
v_e_p = zeros(1,numberPIXI);
v_d_p = zeros(1,numberPIXI);
satellites_p = zeros(1,numberPIXI);
iar_hyp_p = zeros(1,numberPIXI);
iar_ratio_p = zeros(1,numberPIXI);
type_p = zeros(1,numberPIXI);
timeStart = RtkFix.timestamp(1);
%% Extracting PIXI and RTK
j = 1;%PIXI
k = 1;%RTK
fixPixi = 1;
fixRTK = 1;

for i=1:len
    if RtkFix.src_ent(i) == 22
        timestamp_p(j) = RtkFix.timestamp(i);
        src_p(j) = RtkFix.src(i);
        src_ent_p(j) = RtkFix.src_ent(i);
        dst_p(j) =RtkFix.dst(i);
        tow_p(j) = RtkFix.tow(i);
        n_p(j) = RtkFix.n(i);
        e_p(j) = RtkFix.e(i);
        d_p(j) = RtkFix.d(i);
        v_n_p(j) = RtkFix.v_n(i);
        v_e_p(j) =RtkFix.v_e(i);
        v_d_p(j) = RtkFix.v_d(i);
        satellites_p(j) =  RtkFix.satellites(i);
        iar_hyp_p(j) =  RtkFix.iar_hyp(i);
        iar_ratio_p(j) =  RtkFix.iar_ratio(i);
        if strcmp(RtkFix.type(i,1:2),'FI')
            type_p(j) = 3;
%             nf_p(fixPixi) = n_p(j);
%             ef_p(fixPixi) = e_p(j);
%             df_p(fixPixi) = d_p(j);
%             dft_p(fixPixi) = timestamp_p(j);
            fixPixi = fixPixi +1;
        elseif strcmp(RtkFix.type(i,1:2),'FL')
            type_p(j) = 2;
        elseif strcmp(RtkFix.type(i),'O')
            type_p(j) = 1;
        else
            type_p(j) = 0;
        end
        j = j+1;
    else
        timestamp_r(k) = RtkFix.timestamp(i);
        src_r(k) =  RtkFix.src(i);
        src_ent_r(k) =  RtkFix.src_ent(i);
        dst_r(k) = RtkFix.dst(i);
        tow_r(k) =  RtkFix.tow(i);
        n_r(k) =  RtkFix.n(i);
        e_r(k) =  RtkFix.e(i);
        d_r(k) =  RtkFix.d(i);
        v_n_r(k) =  RtkFix.v_n(i);
        v_e_r(k) =  RtkFix.v_e(i);
        v_d_r(k) =  RtkFix.v_d(i);
        satellites_r(k) =  RtkFix.satellites(i);
        iar_hyp_r(k) =  RtkFix.iar_hyp(i);
        iar_ratio_r(k) =  RtkFix.iar_ratio(i);
        if strcmp(RtkFix.type(i,1:2),'FI')
            type_r(k) = 3;
            fixRTK = fixRTK +1;
        elseif strcmp(RtkFix.type(i,1:2),'FL')
            type_r(k) = 2;
        elseif strcmp(RtkFix.type(i), 'O')
            type_r(k) = 1;
        else
            type_r(k) = 0;
        end
        k = k+1;
        
    end
        

end
%% Fixing tow for Piksi, rtklib and PostPro
tow_p = tow_p.*10^-3;
start = 1;
stop = 0;
PostTimeStart = PostPro(:,2) == tow_r(1);
[rowStart,colStart] = find(PostTimeStart);
PostTimeStop = PostPro(:,2) == tow_r(end);
[rowStop,colStop] = find(PostTimeStop);
res = 0.1;
% Since the postprosesing give tow with millisecond, and the output from
% rtkrcv currently does not the closes time need to be found. 
while isempty(rowStart)
    PostTimeStart = abs(PostPro(:,2) - tow_r(1))<=res;
    [rowStart,colStart] = find(PostTimeStart);
    res = res+0.1;
end
res = 0.1;
while isempty(rowStop)
    PostTimeStop = abs(PostPro(:,2) - tow_r(end))<=res;
    [rowStop,colStop] = find(PostTimeStop);
    res = res+0.1;
end
PostTimeStart = PostPro(rowStart(end),2);
PostTimeStop = PostPro(rowStop(end),2);
for i=1:k-1
    if tow_r(i) == tow_r(start)
        stop = stop + 1;
    elseif tow_r(i) ~= tow_r(start)
        
        diff = i-start;
        if diff == 0
            diff = 1;
        end
        for n=0:diff-1
           tow_r(n+start) = tow_r(n+start)+n/diff; 
        end
        start = i;
        stop = i;
    end
end
if tow_r(1)<tow_p(1)
    towTimeStart = tow_r(1);
else
    towTimeStart = tow_p(1);
end
%% Extracting data from PostPro
PostN = PostPro(rowStart(end):rowStop(end),4);
PostE = PostPro(rowStart(end):rowStop(end),3);
PostD = - PostPro(rowStart(end):rowStop(end),5);
PostType = PostPro(rowStart(end):rowStop(end),6);
PostTime = PostPro(rowStart(end):rowStop(end),2);
%% Correcting def off fix,float none from PostPro
for i = 1:length(PostType)
    if PostType(i) == 1 %FIX
        PostType(i) = 3;
    elseif PostType(i) == 2 % FLOAT
        PostType(i) = 2;
    else
        PostType(i) = 0;
    end
end

%% Calculate time difference between timestamp

deltatimetow_r = zeros(1,length(n_r)-1);
deltatimetow_p = zeros(1,length(n_p)-1);
deltatime_r = zeros(1,length(n_r)-1);
deltatime_p = zeros(1,length(n_p)-1);
for i = 1:length(n_r)-1
      deltatime_r(i) = timestamp_r(i+1)-timestamp_r(i);
      deltatimetow_r(i) = tow_r(i+1)-tow_r(i);
end
for i = 1:length(n_p)-1
      deltatime_p(i) = timestamp_p(i+1)-timestamp_p(i);
      deltatimetow_p(i) = tow_p(i+1)-tow_p(i);
end
%% Retrive only fix type solution
% RTK
timestampFix_r = zeros(1,fixRTK-1);
srcFix_r = zeros(1,fixRTK-1);
src_entFix_r = zeros(1,fixRTK-1);
dstFix_r =zeros(1,fixRTK-1);
towFix_r = zeros(1,fixRTK-1);
nFix_r = zeros(1,fixRTK-1);
eFix_r = zeros(1,fixRTK-1);
dFix_r = zeros(1,fixRTK-1);
vFix_n_r = zeros(1,fixRTK-1);
vFix_e_r = zeros(1,fixRTK-1);
vFix_d_r = zeros(1,fixRTK-1);
satellitesFix_r = zeros(1,fixRTK-1);
iar_hypFix_r = zeros(1,fixRTK-1);
iar_ratioFix_r = zeros(1,fixRTK-1);

%Piksi
timestampFix_p = zeros(1,fixPixi-1);
srcFix_p = zeros(1,fixPixi-1);
src_entFix_p = zeros(1,fixPixi-1);
dstFix_p =zeros(1,fixPixi-1);
towFix_p = zeros(1,fixPixi-1);
nFix_p = zeros(1,fixPixi-1);
eFix_p = zeros(1,fixPixi-1);
dFix_p = zeros(1,fixPixi-1);
vFix_n_p = zeros(1,fixPixi-1);
vFix_e_p = zeros(1,fixPixi-1);
vFix_d_p = zeros(1,fixPixi-1);
satellitesFix_p = zeros(1,fixPixi-1);
iar_hypFix_p = zeros(1,fixPixi-1);
iar_ratioFix_p = zeros(1,fixPixi-1);
PiksiFixC = 1;
for i=1:j-1
    if type_p(i)==3
        timestampFix_p(PiksiFixC) = timestamp_p(i);
        srcFix_p(PiksiFixC) = src_p(i);
        src_entFix_p(PiksiFixC) = src_ent_p(i);
        dstFix_p(PiksiFixC) =dst_p(i);
        towFix_p(PiksiFixC) = tow_p(i);
        nFix_p(PiksiFixC) = n_p(i);
        eFix_p(PiksiFixC) = e_p(i);
        dFix_p(PiksiFixC) = d_p(i);
        vFix_n_p(PiksiFixC) = v_n_p(i);
        vFix_e_p(PiksiFixC) = v_e_p(i);
        vFix_d_p(PiksiFixC) = v_d_p(i);
        satellitesFix_p(PiksiFixC) = satellites_p(i);
        iar_hypFix_p(PiksiFixC) = iar_hyp_p(i);
        iar_ratioFix_p(PiksiFixC) = iar_ratio_p(i);
        PiksiFixC = PiksiFixC +1;
    end
end
RTKFixC = 1;
for i=1:k-1
    if type_r(i)==3
        timestampFix_r(RTKFixC) = timestamp_r(i);
        srcFix_r(RTKFixC) = src_r(i);
        src_entFix_r(RTKFixC) = src_ent_r(i);
        dstFix_r(RTKFixC) =dst_r(i);
        towFix_r(RTKFixC) = tow_r(i);
        nFix_r(RTKFixC) = n_r(i);
        eFix_r(RTKFixC) = e_r(i);
        dFix_r(RTKFixC) = d_r(i);
        vFix_n_r(RTKFixC) = v_n_r(i);
        vFix_e_r(RTKFixC) = v_e_r(i);
        vFix_d_r(RTKFixC) = v_d_r(i);
        satellitesFix_r(RTKFixC) = satellites_r(i);
        iar_hypFix_r(RTKFixC) = iar_hyp_r(i);
        iar_ratioFix_r(RTKFixC) = iar_ratio_r(i);
        RTKFixC = RTKFixC +1;
    end
end

%% Creating timeseries
xFix_r = timeseries(nFix_r,towFix_r,'Name','RTKLIB');
yFix_r = timeseries(eFix_r,towFix_r,'Name','RTKLIB');
zFix_r = timeseries(dFix_r,towFix_r,'Name','RTKLIB');
uFix_r = timeseries(vFix_n_r,towFix_r,'Name','RTKLIB');
vFix_r = timeseries(vFix_e_r,towFix_r,'Name','RTKLIB');
wFix_r = timeseries(vFix_d_r,towFix_r,'Name','RTKLIB');

xFix_p = timeseries(nFix_p,towFix_p,'Name','PIKSI');
yFix_p = timeseries(eFix_p,towFix_p,'Name','PIKSI');
zFix_p = timeseries(dFix_p,towFix_p,'Name','PIKSI');
uFix_p = timeseries(vFix_n_p,towFix_p,'Name','PIKSI');
vFix_p = timeseries(vFix_e_p,towFix_p,'Name','PIKSI');
wFix_p = timeseries(vFix_d_p,towFix_p,'Name','PIKSI');

timePostN = timeseries(PostN,PostTime,'Name','PostPro');
timePostE = timeseries(PostE,PostTime,'Name','PostPro');
timePostD = timeseries(PostD,PostTime,'Name','PostPro');
%% Synchornice timeseries

[xpFix_r, timePostN_r] = synchronize(xFix_r,timePostN,'union');
[ypFix_r, timePostE_r] = synchronize(yFix_r,timePostE,'union');
[zpFix_r, timePostD_r] = synchronize(zFix_r,timePostD,'union');
[xpFix_p, timePostN_p] = synchronize(xFix_p,timePostN,'union');
[ypFix_p, timePostE_p] = synchronize(yFix_p,timePostE,'union');
[zpFix_p, timePostD_p] = synchronize(zFix_p,timePostD,'union');


[xFix_r, xFix_p] = synchronize(xFix_r,xFix_p,'union');
[yFix_r, yFix_p] = synchronize(yFix_r,yFix_p,'union');
[zFix_r, zFix_p] = synchronize(zFix_r,zFix_p,'union');
[uFix_r, uFix_p] = synchronize(uFix_r,uFix_p,'union');
[vFix_r, vFix_p] = synchronize(vFix_r,vFix_p,'union');
[wFix_r, wFix_p] = synchronize(wFix_r,wFix_p,'union');

%% Calculate standard deviation of the difference between rtklib and piksi
minLp = length(xpFix_p.data);
minLr = length(xpFix_r.data);
minL = length(xFix_r.Data);
ex = zeros(1,minL);
ey = zeros(1,minL);
ez = zeros(1,minL);
epx_p = zeros(1,minLp);
epy_p = zeros(1,minLp);
epz_p = zeros(1,minLp);
epx_r = zeros(1,minLr);
epy_r = zeros(1,minLr);
epz_r = zeros(1,minLr);
eu = zeros(1,minL);
ev = zeros(1,minL);
ew = zeros(1,minL);

stdx = zeros(1,minL);
stdy = zeros(1,minL);
stdz = zeros(1,minL);
stdxp_p = zeros(1,minLp);
stdyp_p = zeros(1,minLp);
stdzp_p = zeros(1,minLp);
stdxp_r = zeros(1,minLr);
stdyp_r = zeros(1,minLr);
stdzp_r = zeros(1,minLr);
stdu = zeros(1,minL);
stdv = zeros(1,minL);
stdw = zeros(1,minL);

meanx = zeros(1,minL);
meany = zeros(1,minL);
meanz = zeros(1,minL);
meanxp_p = zeros(1,minLp);
meanyp_p = zeros(1,minLp);
meanzp_p = zeros(1,minLp);
meanxp_r = zeros(1,minLr);
meanyp_r = zeros(1,minLr);
meanzp_r = zeros(1,minLr);
meanu = zeros(1,minL);
meanv = zeros(1,minL);
meanw = zeros(1,minL);
for i=1:minL
    ex(i) = xFix_r.Data(i)-xFix_p.Data(i); 
    ey(i) = yFix_r.Data(i)-yFix_p.Data(i);
    ez(i) = zFix_r.Data(i)-wFix_p.Data(i);
    eu(i) = uFix_r.Data(i)-uFix_p.Data(i);
    ev(i) = vFix_r.Data(i)-vFix_p.Data(i);
    ew(i) = wFix_r.Data(i)-wFix_p.Data(i);

    stdx(i) = std(ex(1:i));
    stdy(i) = std(ey(1:i));
    stdz(i) = std(ez(1:i));
    stdu(i) = std(eu(1:i));
    stdv(i) = std(ev(1:i));
    stdw(i) = std(ew(1:i));
    
    meanx(i) = mean(ex(1:i));
    meany(i) = mean(ey(1:i));
    meanz(i) = mean(ez(1:i));
    meanu(i) = mean(eu(1:i));
    meanv(i) = mean(ev(1:i));
    meanw(i) = mean(ew(1:i));

end

for i=1:minLp
    epx_p(i) = xpFix_p.Data(i) - timePostN_p.Data(i);
    epy_p(i) = ypFix_p.Data(i) - timePostE_p.Data(i);
    epz_p(i) = zpFix_p.Data(i) - timePostD_p.Data(i);
    
    stdxp_p(i) = std(epx_p(1:i));
    stdyp_p(i) = std(epy_p(1:i));
    stdzp_p(i) = std(epz_p(1:i));
    
    meanxp_p(i) = mean(epx_p(1:i));
    meanyp_p(i) = mean(epy_p(1:i));
    meanzp_p(i) = mean(epz_p(1:i));
end

for i=1:minLr
    epx_r(i) = xpFix_r.Data(i) - timePostN_r.Data(i);
    epy_r(i) = ypFix_r.Data(i) - timePostE_r.Data(i);
    epz_r(i) = zpFix_r.Data(i) - timePostD_r.Data(i);
    
    stdxp_r(i) = std(epx_r(1:i));
    stdyp_r(i) = std(epy_r(1:i));
    stdzp_r(i) = std(epz_r(1:i));
    
    meanxp_r(i) = mean(epx_r(1:i));
    meanyp_r(i) = mean(epy_r(1:i));
    meanzp_r(i) = mean(epz_r(1:i));
end

%% Plot
TimeEndr = length(timestamp_r);
TimeEndp = length(timestamp_p);
if PIXI
    figure(1);
    plot(eFix_p,nFix_p,'xb');
    hold on;
    plot(eFix_r,nFix_r,'xr');
    plot(PostE,PostN,'--g');
    grid on;
    title('North East position'); 
    xlabel('East [m]'); ylabel('North [m]');
    legend('Piksi','RTKLIB','Post processed');

    % figure(2);
    % plot3(e_p,n_p,d_p);
    % grid on;
    % title('NED Piksi'); 
    % xlabel('East [m]'); ylabel('North [m]'); zlabel('Down [m]');
    % 
    % figure(4);
    % plot3(e_r,n_r,d_r);
    % grid on;
    % title('NED rtklib'); 
    % xlabel('East [m]'); ylabel('North [m]'); zlabel('Down [m]');


    figure(2);
    subplot(2,1,1);
    plot(towFix_p-towFix_p(1),dFix_p,'xb');
    hold on;
    plot(towFix_r-towFix_r(1),dFix_r,'xr');
    plot(PostTime-PostTime(1),PostD,'--g');
    grid on;
    title('Down position');
    legend('Piksi','RTKLIB','Post processed')
    xlabel('Time [s]');
    ylabel('Down [m]');

    subplot(2,1,2)
    plot(tow_p-tow_p(1),type_p,'b');
    hold on;
    plot(tow_r-tow_r(1),type_r,'r');
    plot(PostTime-PostTime(1),PostType,'--g');
    grid on;
    title('Integer ambiguity solution: Fix = 3, Float = 2,Obs = 1, None = 0')
    ylabel('Solution type')
    xlabel('Time [s]');
    legend('Piksi','RTKLIB','Post processed');
    ylim([0 5]);

    figure(3)
    subplot(3,1,1)
    plot(towFix_p-towFix_p(1),vFix_n_p,'--b');
    hold on;
    grid on;
    plot(towFix_r-towFix_r(1),vFix_n_r,'--r')
    legend('Piksi','RTKLIB');
    title('Velocity in North direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,2)
    plot(towFix_p-towFix_p(1),vFix_e_p,'--b');
    hold on;
    grid on;
    plot(towFix_r-towFix_r(1),vFix_e_r,'--r')
    % plot(timestamp_p(1:TimeEnd)-timeStart,ed_p(1:TimeEnd),'g');
    % plot(timestamp_r(1:TimeEnd)-timeStart,ed_r(1:TimeEnd),'c')
    legend('Piksi','RTKLIB');
    title('Velocity in East direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,3)
    plot(timestampFix_p(1:fixPixi-1)-timeStart,vFix_d_p,'--b');
    hold on;
    grid on;
    plot(timestampFix_r(1:fixRTK-1)-timeStart,vFix_d_r,'--r')
    legend('Piksi','RTKLIB');
    title('Velocity in Down direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');
    ylim([-5 5]);

    figure(4)
    plot(timestamp_p(1:TimeEndp)-timeStart,satellites_p(1:TimeEndp),'b');
    hold on;
    plot(timestamp_r(1:TimeEndr)-timeStart,satellites_r(1:TimeEndr),'r');
    grid on;
    title('Visable satellites')
    xlabel('Time [s]');
    ylabel('Number of satellites');
    legend('Piksi','RTKLIB');
    ylim([4 12])
    
    %% Histogram
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
%     
% 	figure(6)
%     subplot(2,1,1);
%     histogram(deltatimetow_p);
%     title('Piksi tow');
%     xlabel('Time [s]');
%     subplot(2,1,2);
%     histogram(deltatime_p);
%     title('Piksi timestamp');
%     xlabel('Time [s]');
%     %% Standard deviation
%     figure(7)
%     plot(xFix_r.Time-towTimeStart,stdx);
%     title('North standard deviation')
%     grid on;
%     xlabel('Time [s]');
%     ylabel('Standard deviation');
%     figure(8)
%     plot(vFix_r.Time-towTimeStart,stdy);
%     title('East standard deviation')
%     grid on;
%     xlabel('Time [s]');
%     ylabel('Standard deviation');
%     figure(9)
%     plot(zFix_r.Time-towTimeStart,stdz);
%     title('Down standard deviation')
%     grid on;
%     xlabel('Time [s]');
%     ylabel('Standard deviation');
%     figure(10)
%     plot(uFix_r.Time-towTimeStart,stdu);
%     title('North velocity standard deviation')
%     grid on;
%     xlabel('Time [s]');
%     ylabel('Standard deviation');
%     figure(11)
%     plot(vFix_r.Time-towTimeStart,stdv);
%     title('East velocity standard deviation')
%     grid on;
%     xlabel('Time [s]');
%     ylabel('Standard deviation');
%     figure(12)
%     plot(wFix_r.Time-towTimeStart,stdw);
%     title('Down velocity standard deviation')
%     grid on;
%     xlabel('Time [s]');
%     ylabel('Standard deviation');
    
    
    
%     figure(13);
%     subplot(2,1,1);
%     plot(towFix_p-towFix_p(1),nFix_p,'b');
%     hold on;
%     plot(towFix_r-towFix_r(1),nFix_r,'r');
%     grid on;
%     title('North Rtklib');
%     legend('Piksi','Rtklib')
%     xlabel('Time [s]');
%     ylabel('Down [m]');
% 
%     subplot(2,1,2)
%     plot(towFix_p-towFix_p(1),eFix_p,'b');
%     hold on;
%     plot(towFix_r-towFix_r(1),eFix_r,'r');
%     grid on;
%     title('East Rtklib');
%     legend('Piksi','Rtklib')
%     xlabel('Time [s]');
%     ylabel('Down [m]');
    
    
    
    %% sync plot
%     figure(14);
%     plot(yFix_p.Data(:),xFix_p.Data(:),'xb');
%     hold on;
%     plot(yFix_r.Data(:),xFix_r.Data(:),'xr');
%     grid on;
%     title('XY'); 
%     xlabel('East [m]'); ylabel('North [m]');
%     legend('Pixi','RtkLib');

    % figure(2);
    % plot3(e_p,n_p,d_p);
    % grid on;
    % title('NED Piksi'); 
    % xlabel('East [m]'); ylabel('North [m]'); zlabel('Down [m]');
    % 
    % figure(4);
    % plot3(e_r,n_r,d_r);
    % grid on;
    % title('NED rtklib'); 
    % xlabel('East [m]'); ylabel('North [m]'); zlabel('Down [m]');


%     figure(15);
%     subplot(2,1,1);
% %     plot(zFix_p.Time-timeStart,zFix_p.Data(:),'xb');
%     plot(zpFix_p.Time-towTimeStart,zpFix_p.Data(:),'xb');
%     hold on;
% %     plot(zFix_r.Time-timeStart,zFix_r.Data(:),'xr');
%     plot(timePostD_p.Time -towTimeStart,timePostD_p.Data(:),'xr');
%     grid on;
%     title('Down Rtklib');
%     legend('Pixi','Rtklib')
%     xlabel('Time [s]');
%     ylabel('Down [m]');

%     subplot(2,1,2)
%     plot(timestamp_p(1:TimeEndp)-timeStart,type_p(1:TimeEndp),'b');
%     hold on;
%     plot(timestamp_r(1:TimeEndr)-timeStart,type_r(1:TimeEndr),'r');
%     grid on;
%     title('Ambiguity solution')
%     ylabel('Solution type')
%     xlabel('Time [s]');
%     legend('Pixi','Rtklib');
%     ylim([0 5]);

%     figure(16)
%     subplot(3,1,1)
%     plot(uFix_p.Time-timeStart,uFix_p.Data(:),'xb');
%     hold on;
%     grid on;
%     plot(uFix_r.Time-timeStart,uFix_r.Data(:),'xr')
%     legend('Pixi','Rtklib');
%     title('Velocity in North direction');
%     ylabel('Velocity [m/s]')
%     xlabel('Time [s]');
% 
%     subplot(3,1,2)
%     plot(vFix_p.Time-timeStart,vFix_p.Data(:),'xb');
%     hold on;
%     grid on;
%     plot(vFix_r.Time-timeStart,vFix_r.Data(:),'xr')
%     % plot(timestamp_p(1:TimeEnd)-timeStart,ed_p(1:TimeEnd),'g');
%     % plot(timestamp_r(1:TimeEnd)-timeStart,ed_r(1:TimeEnd),'c')
%     legend('Pixi','Rtklib');
%     title('Velocity in East direction');
%     ylabel('Velocity [m/s]')
%     xlabel('Time [s]');
% 
%     subplot(3,1,3)
%     plot(wFix_p.Time-timeStart,wFix_p.Data(:),'xb');
%     hold on;
%     grid on;
%     plot(wFix_r.Time-timeStart,wFix_r.Data(:),'xr')
%     legend('Pixi','Rtklib');
%     title('Velocity in Down direction');
%     ylabel('Velocity [m/s]')
%     xlabel('Time [s]');
    
    %% Mean plot
%     figure(17);
%     plot(xFix_r.Time-towTimeStart,meanx);
%     grid on;
%     title('The mean difference in North estimation');
%     ylabel('Mean difference [m]');
%     xlabel('Time [s]');
%     figure(18);
%     plot(yFix_r.Time-towTimeStart,meany);
%     grid on;
%     title('The mean difference in East estimation');
%     ylabel('Mean difference [m]');
%     xlabel('Time [s]');
%     figure(19);
%     plot(zFix_r.Time-towTimeStart,meanz);
%     grid on;
%     title('The mean difference in Down estimation');
%     ylabel('Mean difference [m]');
%     xlabel('Time [s]');
%     figure(20);
%     plot(uFix_r.Time-towTimeStart,meanu);
%     grid on;
%     title('The mean difference in North velocity estimation');
%     ylabel('Mean difference [m/s]');
%     xlabel('Time [s]');
%     figure(21);
%     plot(vFix_r.Time-towTimeStart,meanv);
%     grid on;
%     title('The mean difference in East velocity estimation');
%     ylabel('Mean difference [m/s]');
%     xlabel('Time [s]');
%     figure(22);
%     plot(wFix_r.Time-towTimeStart,meanw);
%     grid on;
%     title('The mean difference in Down velocity estimation');
%     ylabel('Mean difference [m/s]');
%     xlabel('Time [s]');
    
    %% Plot were the post proseccing is used
    figure(23)
    subplot(3,1,1)
    plot(xpFix_r.Time - towTimeStart,epx_r);
    grid on;
    title('North error, RTKLIB');
    ylabel('Distance [m]');
    xlabel('Time [s]');
    subplot(3,1,2)
    plot(ypFix_r.Time - towTimeStart,epy_r);
    grid on;
    title('East error, RTKLIB');
    ylabel('Distance [m]');
    xlabel('Time [s]');
    subplot(3,1,3)
    plot(zpFix_r.Time - towTimeStart,epz_r);
    grid on;
    title('Down error, RTKLIB');
    ylabel('Distance [m]');
    xlabel('Time [s]');
    figure(24)
    subplot(3,1,1);
    plot(xpFix_r.Time-towTimeStart,stdxp_r);
    grid on;
    title('Cumulative standard deviation in North position, RTKLIB');
    ylabel('Distance [m]');
    xlabel('Time [s]');
    subplot(3,1,2);
    plot(ypFix_r.Time-towTimeStart,stdyp_r);
    grid on;
    title('Cumulative standard deviation in East position, RTKLIB');
    ylabel('Distance [m]');
    xlabel('Time [s]');
    subplot(3,1,3);
    plot(zpFix_r.Time-towTimeStart,stdzp_r);
    grid on;
    title('Cumulative standard deviation in Down position, RTKLIB');
    ylabel('Distance [m]');
    xlabel('Time [s]');
%     figure(25)
%     subplot(3,1,1);
%     plot(xpFix_r.Time-towTimeStart,meanxp_r);
%     grid on;
%     title('Mean value Rtklib and Post');
%     ylabel('Distance [m]');
%     xlabel('Time [s]');
%     subplot(3,1,2);
%     plot(ypFix_r.Time-towTimeStart,meanyp_r);
%     grid on;
%     title('Mean value Rtklib and Post');
%     ylabel('Distance [m]');
%     xlabel('Time [s]');
%     subplot(3,1,3);
%     plot(zpFix_r.Time-towTimeStart,meanzp_r);
%     grid on;
%     title('Mean value Rtklib and Post');
%     ylabel('Distance [m]');
%     xlabel('Time [s]');
    
    figure(26);
    subplot(3,1,1)
    plot(xpFix_p.Time - towTimeStart,epx_p);
    grid on;
    title('North error, Piksi');
    ylabel('Distance [m]');
    xlabel('Time [s]');
    subplot(3,1,2)
    plot(ypFix_p.Time - towTimeStart,epy_p);
    grid on;
    title('East error, Piksi');
    ylabel('Distance [m]');
    xlabel('Time [s]');
    subplot(3,1,3)
    plot(zpFix_p.Time - towTimeStart,epz_p);
    grid on;
    title('Down error, Piksi');
    ylabel('Distance [m]');
    xlabel('Time [s]');
    figure(27)
    subplot(3,1,1);
    plot(xpFix_p.Time-towTimeStart,stdxp_p);
    grid on;
    title('Cumulative standard deviation in North position, Piksi');
    ylabel('Distance [m]');
    xlabel('Time [s]');
    subplot(3,1,2);
    plot(ypFix_p.Time-towTimeStart,stdyp_p);
    grid on;
    title('Cumulative standard deviation in East position, Piksi');
    ylabel('Distance [m]');
    xlabel('Time [s]');
    subplot(3,1,3);
    plot(zpFix_p.Time-towTimeStart,stdzp_p);
    grid on;
    title('Cumulative standard deviation in Down position, Piksi');
    ylabel('Distance [m]');
    xlabel('Time [s]');
%     figure(28)
%     subplot(3,1,1);
%     plot(xpFix_p.Time-towTimeStart,meanxp_p);
%     grid on;
%     title('Cumulative mean value Piksi and Post');
%     ylabel('Distance [m]');
%     xlabel('Time [s]');
%     subplot(3,1,2);
%     plot(ypFix_p.Time-towTimeStart,meanyp_p);
%     grid on;
%     title('Cumulative mean value Piksi and Post');
%     ylabel('Distance [m]');
%     xlabel('Time [s]');
%     subplot(3,1,3);
%     plot(zpFix_p.Time-towTimeStart,meanzp_p);
%     grid on;
%     title('Cumulative mean value Piksi and Post');
%     ylabel('Distance [m]');
%     xlabel('Time [s]');
    figure(29);
    plot(towFix_p-towFix_p(1),dFix_p,'b');
    hold on;
    plot(towFix_r-towFix_r(1),dFix_r,'r');
    plot(PostTime-PostTime(1),PostD,'--g');
    grid on;
    title('Down position');
    legend('Piksi','RTKLIB','Post processed')
    xlabel('Time [s]');
    ylabel('Position [m]');  
else
    figure(1);
    plot(e_r,n_r,'r');
    grid on;
    title('XY'); 
    xlabel('East [m]'); ylabel('North [m]');

    % figure(2);
    % plot3(e_p,n_p,d_p);
    % grid on;
    % title('NED Piksi'); 
    % xlabel('East [m]'); ylabel('North [m]'); zlabel('Down [m]');
    % 
    % figure(4);
    % plot3(e_r,n_r,d_r);
    % grid on;
    % title('NED rtklib'); 
    % xlabel('East [m]'); ylabel('North [m]'); zlabel('Down [m]');


    figure(2);
    subplot(2,1,1);
    plot(timestamp_r(1:TimeEndr)-timeStart,d_r(1:TimeEndr),'r');
    grid on;
    title('Down Rtklib');
    xlabel('Time [s]');
    ylabel('Down [m]');

    subplot(2,1,2)
    plot(timestamp_r(1:TimeEndr)-timeStart,type_r(1:TimeEndr),'r');
    grid on;
    title('Abiguity solution')
    ylabel('Solution type')
    xlabel('Time [s]');
    ylim([0 5]);

    figure(3)
    subplot(3,1,1)
    grid on;
    plot(timestamp_r(1:TimeEndr)-timeStart,v_n_r(1:TimeEndr),'r')
    title('Velocity in North direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,2)
    grid on;
    plot(timestamp_r(1:TimeEndr)-timeStart,v_e_r(1:TimeEndr),'r')
    title('Velocity in East direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,3)
    grid on;
    plot(timestamp_r(1:TimeEndr)-timeStart,v_d_r(1:TimeEndr),'r')
    title('Velocity in Down direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    figure(4)
    plot(timestamp_r(1:TimeEndr)-timeStart,satellites_r(1:TimeEndr),'r');
    grid on;
    title('Visable satellites')
    xlabel('Time [s]');
    ylabel('Number of satellites');
    ylim([4 12])
end

