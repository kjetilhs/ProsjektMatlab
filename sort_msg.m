%piksi 22
%rtklib 23

load('Data.mat')

len = length(RtkFix.src_ent);

timestamp_r = [];
src_r = [];
src_ent_r = [];
dst_r =[];
tow_r = [];
n_r = [];
e_r = [];
d_r = [];
v_n_r = [];
v_e_r = [];
v_d_r = [];
satellites_r = [];
iar_hyp_r = [];
iar_ratio_r = [];
type_r = [];

timestamp_p = [];
src_p = [];
src_ent_p = [];
dst_p =[];
tow_p = [];
n_p = [];
e_p = [];
d_p = [];
v_n_p = [];
v_e_p = [];
v_d_p = [];
satellites_p = [];
iar_hyp_p = [];
iar_ratio_p = [];
type_p = [];


for i=1:len
    if RtkFix.src_ent(i) == 22
        timestamp_p = [timestamp_p  RtkFix.timestamp(i)];
        src_p = [src_p RtkFix.src(i)];
        src_ent_p = [src_ent_p RtkFix.src_ent(i)];
        dst_p =[dst_p RtkFix.dst(i)];
        tow_p = [tow_p RtkFix.tow(i)];
        n_p = [n_p RtkFix.n(i)];
        e_p = [e_p RtkFix.e(i)];
        d_p = [d_p RtkFix.d(i)];
        v_n_p = [v_n_p RtkFix.v_n(i)];
        v_e_p = [v_e_p RtkFix.v_e(i)];
        v_d_p = [v_d_p RtkFix.v_d(i)];
        satellites_p = [satellites_p RtkFix.satellites(i)];
        iar_hyp_p = [iar_hyp_p RtkFix.iar_hyp(i)];
        iar_ratio_p = [iar_ratio_p RtkFix.iar_ratio(i)];
        type_p = [type_p RtkFix.type(i)];
    
    else
        timestamp_r = [timestamp_r  RtkFix.timestamp(i)];
        src_r = [src_r RtkFix.src(i)];
        src_ent_r = [src_ent_r RtkFix.src_ent(i)];
        dst_r =[dst_r RtkFix.dst(i)];
        tow_r = [tow_r RtkFix.tow(i)];
        n_r = [n_r RtkFix.n(i)];
        e_r = [e_r RtkFix.e(i)];
        d_r = [d_r RtkFix.d(i)];
        v_n_r = [v_n_r RtkFix.v_n(i)];
        v_e_r = [v_e_r RtkFix.v_e(i)];
        v_d_r = [v_d_r RtkFix.v_d(i)];
        satellites_r = [satellites_r RtkFix.satellites(i)];
        iar_hyp_r = [iar_hyp_r RtkFix.iar_hyp(i)];
        iar_ratio_r = [iar_ratio_r RtkFix.iar_ratio(i)];
        type_r = [type_r RtkFix.type(i)];
        
        
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
title('ENU Piksi'); 
xlabel('East [m]'); ylabel('North [m]'); zlabel('Up [m]');

figure(4);
plot3(e_r,n_r,d_r);
grid on;
title('ENU rtklib'); 
xlabel('East [m]'); ylabel('North [m]'); zlabel('Up [m]');


figure(5);
plot(timestamp_p,d_p);
grid on;
title('Down Piksi'); 
xlabel('Time [s]'); ylabel('Down [m]');

figure(6);
plot(timestamp_r,d_r);
grid on;
title('Down Rtklib'); 
xlabel('Time [s]'); ylabel('Down [m]');