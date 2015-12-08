PostPro = load('log_folder/PostPro/rtklib_ubxstream_x8_log201512081128.pos');
flightdata = load('log_folder/rtklibLog/rtklib_output201512081128.pos');
len_msg = length(flightdata(:,1)); 

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
vd = flightdata(tStart:tStop,18);

tp = flightdata(tStart:tStop,2)-flightdata(tStart,2);
ep = flightdata(tStart:tStop,3);
np = flightdata(tStart:tStop,4);
dp = -flightdata(tStart:tStop,5);
qualityp = flightdata(tStart:tStop,6);
num_satp = flightdata(tStart:tStop,7);
sdxp = flightdata(tStart:tStop,8);
sdyp = flightdata(tStart:tStop,9);
sdzp = flightdata(tStart:tStop,10);
sdxyp = flightdata(tStart:tStop,11);
sdyzp = flightdata(tStart:tStop,12);
sdzxp = flightdata(tStart:tStop,13);
agep = flightdata(tStart:tStop,14);
ratiop = flightdata(tStart:tStop,15);


tStart = 1542;
tStop = 4958;
figure(1);
    plot(eFix_p,nFix_p,'xb');
    hold on;
    plot(eFix_r,nFix_r,'xr');
    plot(PostE,PostN,'xg');
    grid on;
    title('XY'); 
    xlabel('East [m]'); ylabel('North [m]');
    legend('Pixi','RtkLib','post processed rtklib');

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
    plot(PostTime-PostTime(1),PostD,'xg');
    grid on;
    title('Down');
    legend('Pixi','Rtklib','post processed rtklib')
    xlabel('Time [s]');
    ylabel('Down [m]');

    subplot(2,1,2)
    plot(tow_p-tow_p(1),type_p,'b');
    hold on;
    plot(tow_r-tow_r(1),type_r,'r');
    plot(PostTime-PostTime(1),PostType,'--g');
    grid on;
    title('Ambiguity solution')
    ylabel('Solution type')
    xlabel('Time [s]');
    legend('Pixi','Rtklib','post processed rtklib');
    ylim([0 5]);

    figure(3)
    subplot(3,1,1)
    plot(timestampFix_p(1:fixPixi-1)-timeStart,vFix_n_p(1:fixPixi-1),'xb');
    hold on;
    grid on;
    plot(timestampFix_r(1:fixRTK-1)-timeStart,vFix_n_r(1:fixRTK-1),'xr')
    legend('Pixi','Rtklib');
    title('Velocity in North direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,2)
    plot(timestampFix_p(1:fixPixi-1)-timeStart,vFix_e_p(1:fixPixi-1),'xb');
    hold on;
    grid on;
    plot(timestampFix_r(1:fixRTK-1)-timeStart,v_e_r(1:fixRTK-1),'xr')
    % plot(timestamp_p(1:TimeEnd)-timeStart,ed_p(1:TimeEnd),'g');
    % plot(timestamp_r(1:TimeEnd)-timeStart,ed_r(1:TimeEnd),'c')
    legend('Pixi','Rtklib');
    title('Velocity in East direction');
    ylabel('Velocity [m/s]')
    xlabel('Time [s]');

    subplot(3,1,3)
    plot(timestampFix_p(1:fixPixi-1)-timeStart,vFix_d_p(1:fixPixi-1),'xb');
    hold on;
    grid on;
    plot(timestampFix_r(1:fixRTK-1)-timeStart,vFix_d_r(1:fixRTK-1),'xr')
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
    
    %% Histogram
    figure(5)
    subplot(2,1,1);
    histogram(deltatimetow_r)
    title('RTKLIB tow');
    xlabel('Time [s]');
    
    subplot(2,1,2);
    histogram(deltatime_r)
    title('RTKLIB timestamp');
    xlabel('Time [s]');