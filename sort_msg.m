%piksi 22
%rtklib 23

load('log_folder/neptusLog/Springeoekt1/Data.mat')

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
        type_p(j) =  RtkFix.type(i);
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
        type_r(k) =  RtkFix.type(i);
        k = k+1;
        
    end
        
        

end

figure(1);
plot(e_p,n_p);
grid on;
title('Piksi'); 
xlabel('East [m]'); ylabel('North [m]');

figure(2);
plot(e_r,n_r);
grid on;
title('rtklib'); 
xlabel('East [m]'); ylabel('North [m]');

figure(3);
plot3(e_p,n_p,d_p);
grid on;
title('NED Piksi'); 
xlabel('East [m]'); ylabel('North [m]'); zlabel('Down [m]');

figure(4);
plot3(e_r,n_r,d_r);
grid on;
title('NED rtklib'); 
xlabel('East [m]'); ylabel('North [m]'); zlabel('Down [m]');


figure(5);
plot(timestamp_p(1:190)-timeStart,d_p(1:190));
grid on;
title('Down Piksi'); 
xlabel('Time [s]'); ylabel('Down [m]');

figure(6);
plot(timestamp_r(1:190)-timeStart,d_r(1:190));
grid on;
title('Down Rtklib'); 
xlabel('Time [s]'); ylabel('Down [m]');