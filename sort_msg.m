%piksi 22
%rtklib 23
close all;
clear;
load('log_folder/neptusLog/Vandring4/Data.mat')
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
%% Calculate time difference between timestamp

deltatime_r = zeros(1,length(n_r)-1);
deltatime_p = zeros(1,length(n_p)-1);

for i = 1:length(n_r)-1
    deltatime_r(i) = timestamp_r(i+1)-timestamp_r(i);
end
for i = 1:length(n_p)-1
    deltatime_p(i) = timestamp_p(i+1)-timestamp_p(i);
end
%% Retrive only fix type solution
% RTK
timestampFix_r = zeros(1,fixRTK);
srcFix_r = zeros(1,fixRTK);
src_entFix_r = zeros(1,fixRTK);
dstFix_r =zeros(1,fixRTK);
towFix_r = zeros(1,fixRTK);
nFix_r = zeros(1,fixRTK);
eFix_r = zeros(1,fixRTK);
dFix_r = zeros(1,fixRTK);
vFix_n_r = zeros(1,fixRTK);
vFix_e_r = zeros(1,fixRTK);
vFix_d_r = zeros(1,fixRTK);
satellitesFix_r = zeros(1,fixRTK);
iar_hypFix_r = zeros(1,fixRTK);
iar_ratioFix_r = zeros(1,fixRTK);

%Piksi
timestampFix_p = zeros(1,fixPixi);
srcFix_p = zeros(1,fixPixi);
src_entFix_p = zeros(1,fixPixi);
dstFix_p =zeros(1,fixPixi);
towFix_p = zeros(1,fixPixi);
nFix_p = zeros(1,fixPixi);
eFix_p = zeros(1,fixPixi);
dFix_p = zeros(1,fixPixi);
vFix_n_p = zeros(1,fixPixi);
vFix_e_p = zeros(1,fixPixi);
vFix_d_p = zeros(1,fixPixi);
satellitesFix_p = zeros(1,fixPixi);
iar_hypFix_p = zeros(1,fixPixi);
iar_ratioFix_p = zeros(1,fixPixi);
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

%% Calculate standard deviation of the difference between rtklib and piksi
minL = min([fixPixi fixRTK]);
ex = zeros(1,minL);
ey = zeros(1,minL);
ez = zeros(1,minL);
eu = zeros(1,minL);
ev = zeros(1,minL);
ew = zeros(1,minL);

stdx = zeros(1,minL);
stdy = zeros(1,minL);
stdz = zeros(1,minL);
stdu = zeros(1,minL);
stdv = zeros(1,minL);
stdw = zeros(1,minL);

for i=1:minL
    ex(i) = nFix_r(i)-nFix_p(i); 
    ey(i) = eFix_r(i)-eFix_p(i);
    ez(i) = dFix_r(i)-dFix_p(i);
    eu(i) = vFix_n_r(i)-vFix_n_p(i);
    ev(i) = vFix_e_r(i)-vFix_e_p(i);
    ew(i) = vFix_d_r(i)-vFix_d_p(i);

    stdx(i) = std(ex(1:i));
    stdy(i) = std(ey(1:i));
    stdz(i) = std(ez(1:i));
    stdu(i) = std(eu(1:i));
    stdv(i) = std(ev(1:i));
    stdw(i) = std(ew(1:i));
end

%% Plot
TimeEndr = length(timestamp_r);
TimeEndp = length(timestamp_p);
if PIXI
    figure(1);
    plot(e_p,n_p,'xb');
    hold on;
    plot(e_r,n_r,'xr');
    grid on;
    title('XY'); 
    xlabel('East [m]'); ylabel('North [m]');
    legend('Pixi','RtkLib');

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
    plot(timestamp_p(1:TimeEndp)-timeStart,d_p(1:TimeEndp),'xb');
%     plot(dft_p-timeStart,df_p,'b');
    % grid on;
    % title('Down Piksi'); 
    % xlabel('Time [s]'); ylabel('Down [m]');

    hold on;
    plot(timestamp_r(1:TimeEndr)-timeStart,d_r(1:TimeEndr),'xr');
    grid on;
    title('Down Rtklib');
    legend('Pixi','Rtklib')
    xlabel('Time [s]');
    ylabel('Down [m]');

    subplot(2,1,2)
    plot(timestamp_p(1:TimeEndp)-timeStart,type_p(1:TimeEndp),'b');
    hold on;
    plot(timestamp_r(1:TimeEndr)-timeStart,type_r(1:TimeEndr),'r');
    grid on;
    title('Ambiguity solution')
    ylabel('Solution type')
    xlabel('Time [s]');
    legend('Pixi','Rtklib');
    ylim([0 5]);

    figure(3)
    subplot(3,1,1)
    plot(timestamp_p(1:TimeEndp)-timeStart,v_n_p(1:TimeEndp),'xb');
    hold on;
    grid on;
    plot(timestamp_r(1:TimeEndr)-timeStart,v_n_r(1:TimeEndr),'xr')
    legend('Pixi','Rtklib');
    title('Velocity in North direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,2)
    plot(timestamp_p(1:TimeEndp)-timeStart,v_e_p(1:TimeEndp),'xb');
    hold on;
    grid on;
    plot(timestamp_r(1:TimeEndr)-timeStart,v_e_r(1:TimeEndr),'xr')
    % plot(timestamp_p(1:TimeEnd)-timeStart,ed_p(1:TimeEnd),'g');
    % plot(timestamp_r(1:TimeEnd)-timeStart,ed_r(1:TimeEnd),'c')
    legend('Pixi','Rtklib');
    title('Velocity in East direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,3)
    plot(timestamp_p(1:TimeEndp)-timeStart,v_d_p(1:TimeEndp),'xb');
    hold on;
    grid on;
    plot(timestamp_r(1:TimeEndr)-timeStart,v_d_r(1:TimeEndr),'xr')
    legend('Pixi','Rtklib');
    title('Velocity in Down direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    figure(4)
    plot(timestamp_p(1:TimeEndp)-timeStart,satellites_p(1:TimeEndp),'b');
    hold on;
    plot(timestamp_r(1:TimeEndr)-timeStart,satellites_r(1:TimeEndr),'r');
    grid on;
    title('Visable satellites')
    xlabel('Time [s]');
    ylabel('Number of satellites');
    legend('Pixi','Rtklib');
    ylim([4 12])
    
    
    figure(5) 
    histogram(deltatime_r)
    title('RTKLIB');
    xlabel('Time [s]');
	figure(6)
    histogram(deltatime_p);
    title('Piksi');
    figure(7)
    plot(stdx);
    figure(8)
    plot(stdy);
    figure(9)
    plot(stdz);
    figure(10)
    plot(stdu);
    figure(11)
    plot(stdv);
    figure(12)
    plot(stdw);
    
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

