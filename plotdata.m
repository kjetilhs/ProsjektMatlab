flightdata = load('rtklib_output2000010101.pos');

len_msg = length(flightdata(:,1)); 

t = flightdata(:,2)-flightdata(1,2);
e = flightdata(:,3);
n = flightdata(:,4);
u = flightdata(:,5);
quality = flightdata(:,6);
num_sat = flightdata(:,7);
sdx = flightdata(:,8);
sdy = flightdata(:,9);
sdz = flightdata(:,10);
sdxy = flightdata(:,11);
sdyz = flightdata(:,12);
sdzx = flightdata(:,13);
age = flightdata(:,14);
ratio = flightdata(:,15);
%{
v_e = flightdata(:,16);
v_n = flightdata(:,17);
v_u = flightdata(:,18);
%}

figure(1);
plot(e,n);
grid on;
title('ENU plot of estimated postion'); 
xlabel('East [m]'); ylabel('North [m]');

figure(2);
plot3(e,n,u);
grid on;
title('ENU plot of estimated postion'); 
xlabel('East [m]'); ylabel('North [m]'); zlabel('Up [m]');

figure(3);
plot(t,u);
grid on;
title('Up plot of estimated position'); 
xlabel('Time [s]'); ylabel('Up [m]');

figure(4);
plot3(n,e,-u);
grid on;
title('NED plot of estimated postion'); 
xlabel('North [m]'); ylabel('East [m]'); zlabel('Down [m]');

figure(5);
plot(t,num_sat);
grid on;
axis([0 t(end) 5 15]);
title('Number of satellites during flight'); 
xlabel('Time [s]'); ylabel('Number of satellites');

figure(6);
plot(t,quality)
axis([0 t(end) 0 5]);
title('Quality'); 
xlabel('Time [s]'); ylabel('Quality of the position data');

figure(7);
subplot(3,1,1);
plot(t,sdx);
grid on
title('Standard Deviation'); ylabel('x');
subplot(3,1,2);
plot(t,sdy);
grid on
ylabel('y');
subplot(3,1,3);
plot(t,sdz);
grid on
ylabel('z');

figure(8);
subplot(3,1,1);
plot(t,sdxy);
grid on
title('Standard Deviation'); ylabel('xy');
subplot(3,1,2);
plot(t,sdyz);
grid on
ylabel('yz');
subplot(3,1,3);
plot(t,sdzx);
grid on
ylabel('zx');










