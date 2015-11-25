%piksi 22
%rtklib 23
close all;
clear;
load('log_folder/neptusLog/Vandring3/Data.mat')
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
            type_p(k) = 3;
        elseif strcmp(RtkFix.type(i,1:2),'FL')
            type_p(k) = 2;
        elseif strcmp(RtkFix.type(i),'O')
            type_p(k) = 1;
        else
            type_p(k) = 0;
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
%% Calculate velocity from position
nd_r = zeros(1,length(n_r)-1);
ed_r = zeros(1,length(n_r)-1);
dd_r = zeros(1,length(n_r)-1);
nd_p = zeros(1,length(n_p)-1);
ed_p = zeros(1,length(n_p)-1);
dd_p = zeros(1,length(n_p)-1);
for i = 1:length(n_r)-1
    nd_r(i) = (n_r(i+1)-n_r(i))/(timestamp_r(i+1)-timestamp_r(i));
    ed_r(i) = (e_r(i+1)-e_r(i))/(timestamp_r(i+1)-timestamp_r(i));
    dd_r(i) = (d_r(i+1)-d_r(i))/(timestamp_r(i+1)-timestamp_r(i));
end
for i = 1:length(n_p)-1
    nd_p(i) = (n_p(i+1)-n_p(i))/(timestamp_p(i+1)-timestamp_p(i));
    ed_p(i) = (e_p(i+1)-e_p(i))/(timestamp_p(i+1)-timestamp_p(i));
    dd_p(i) = (d_p(i+1)-d_p(i))/(timestamp_p(i+1)-timestamp_p(i));
end
%% Plot
TimeEndr = length(timestamp_r);
TimeEndp = length(timestamp_p);
if PIXI
    figure(1);
    plot(e_p,n_p);
    hold on;
    plot(e_r,n_r,'r');
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
    plot(timestamp_p(1:TimeEndp)-timeStart,d_p(1:TimeEndp),'b');
    % grid on;
    % title('Down Piksi'); 
    % xlabel('Time [s]'); ylabel('Down [m]');

    hold on;
    plot(timestamp_r(1:TimeEndr)-timeStart,d_r(1:TimeEndr),'r');
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
    title('Abiguity solution')
    ylabel('Solution type')
    xlabel('Time [s]');
    legend('Pixi','Rtklib');
    ylim([0 5]);

    figure(3)
    subplot(3,1,1)
    plot(timestamp_p(1:TimeEndp)-timeStart,v_n_p(1:TimeEndp));
    hold on;
    grid on;
    plot(timestamp_r(1:TimeEndr)-timeStart,v_n_r(1:TimeEndr),'r')
    legend('Pixi','Rtklib');
    title('Velocity in North direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,2)
    plot(timestamp_p(1:TimeEndp)-timeStart,v_e_p(1:TimeEndp));
    hold on;
    grid on;
    plot(timestamp_r(1:TimeEndr)-timeStart,v_e_r(1:TimeEndr),'r')
    % plot(timestamp_p(1:TimeEnd)-timeStart,ed_p(1:TimeEnd),'g');
    % plot(timestamp_r(1:TimeEnd)-timeStart,ed_r(1:TimeEnd),'c')
    legend('Pixi','Rtklib');
    title('Velocity in East direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,3)
    plot(timestamp_p(1:TimeEndp)-timeStart,v_d_p(1:TimeEndp));
    hold on;
    grid on;
    plot(timestamp_r(1:TimeEndr)-timeStart,v_d_r(1:TimeEndr),'r')
    legend('Pixi','Rtklib');
    title('Velocity in Down direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    figure(4)
    plot(timestamp_p(1:TimeEndp)-timeStart,satellites_p(1:TimeEndp));
    hold on;
    plot(timestamp_r(1:TimeEndr)-timeStart,satellites_r(1:TimeEndr),'r');
    grid on;
    title('Visable satellites')
    xlabel('Time [s]');
    ylabel('Number of satellites');
    legend('Pixi','Rtklib');
    ylim([4 12])
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