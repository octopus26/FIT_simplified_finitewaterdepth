%clear all;
close all;
%% Parameters timed by scale ratio of 30
Buoy_R=3; %[m]
Buoy_Draft_o= 43.2; %[m]

% Buoy_R= 0.1; %[m]
% Buoy_Draft_o= 1.44; %[m]
gravity = 9.80665; 
rho=1027;
Cd = 0.7;

z = linspace(0,(-1)*Buoy_Draft_o,500);
dz = z(1)-z(2);
x = 0;
%% Read input wave elevation from experiment
%for i = 1:15
[domega,omega,amp_s,phase,Fx,My,time,wave,dt] = input_exp_full;

amp_s = amp_s*1/pi*domega;
omega_s = omega;
time_record = time;
fileID4 = fopen('time.txt','at');
fprintf(fileID4,'\n%e',time');
% A(omega) & k(omega)
%k_w_s =omega_s.^2./gravity;  %deep water
dispersion = @(k_w_s) omega_s.^2 - k_w_s*gravity.*tanh(k_w_s*Buoy_Draft_o);
k_w_s = fsolve(dispersion,ones(1,length(omega_s)));

fileID1 = fopen('amp.txt','at');
fprintf(fileID1,'\n%e',amp_s');
fileID2 = fopen('omega.txt','at');
fprintf(fileID2,'\n%e',omega_s');
%Incident wave potential 
% phase is [-pi,pi]
fileID3 = fopen('phase.txt','at');
fprintf(fileID3,'\n%e',phase');
% Phio = zeros(length(time_record),length(z));
% for count_z = 1:length(z)
%     for count_t = 1:length(time_record)
%         Phio(count_t,count_z) = real(sum(1i*gravity*amp_s./omega_s.*exp(k_w_s.*z(count_z)-1i*k_w_s.*x + 1i*omega_s.*time_record(count_t)+ 1i*phase)));
%     end 
% end 
%%
%Incident wave elevation 
zetao = zeros(1,length(time_record));
for count_t = 1:length(time_record)
zetao(count_t) = real(sum(amp_s.*exp(-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
end

%%
%Incident wave velocity
u1 = zeros(length(time_record),length(z));
u1dot = zeros(length(time_record),length(z));
u1x = zeros(length(time_record),length(z));
u1z = zeros(length(time_record),length(z));
u3 = zeros(length(time_record),length(z));
for count_z = 1:length(z)
    for count_t = 1:length(time_record)
        u1(count_t,count_z) = real(sum((-1i)*k_w_s.*1i*gravity.*amp_s./omega_s.*exp(k_w_s.*z(count_z)-1i*k_w_s.*x + 1i*omega_s.*time_record(count_t)+ 1i*phase)));
        u1dot(count_t,count_z) = real(sum((1i*omega_s).*((-1i)*k_w_s).*1i*gravity.*amp_s./omega_s.*exp(k_w_s*z(count_z)-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
        u1x(count_t,count_z) = real(sum(((-1i)*k_w_s).^2.*1i*gravity.*amp_s./omega_s.*exp(k_w_s*z(count_z)-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
        u1z(count_t,count_z) = real(sum(k_w_s.*((-1i)*k_w_s).*1i*gravity.*amp_s./omega_s.*exp(k_w_s*z(count_z)-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
        u3(count_t,count_z) = real(sum(k_w_s.*1i*gravity.*amp_s./omega_s.*exp(k_w_s*z(count_z)-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
       
    end 
end 

% %Derivatives of incident wave elevation
zetaot = zeros(1,length(count_t));
zetaox = zeros(1,length(count_t));
zetaotx = zeros(1,length(count_t));
for count_t = 1:length(time_record)
zetaot(count_t) = real(sum(1i*omega_s.*amp_s.*exp(-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
zetaox(count_t) = real(sum(-1i*k_w_s.*amp_s.*exp(-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
zetaotx(count_t) = real(sum(amp_s.*omega_s.*k_w_s.*exp(-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
end

%% Morison Force
load('dphi2dxdt');
u2dot = dphi2dxdt;
F_Morison = 2*rho*pi*Buoy_R^2*(trapz(u1dot,2)')*dz;  %Classic Morison 
%F_Morison_2 = 2*rho*pi*Buoy_R^2*(trapz(u2dot,2)')*dz;
F_body_nonlinear = 2*rho*pi*Buoy_R^2*trapz(u1.*u1x+u3.*u1z,2)'*dz;
F_wave_ele_nonl = 2*rho*pi*Buoy_R^2*zetao.*(u1dot(:,1)+u1(:,1).*u1x(:,1)+u3(:,1).*u1z(:,1))';

F_crest = rho*pi*Buoy_R^2*u1(:,1)'.*(zetaot-u1(:,1)'.*zetaox);
F_second = rho*pi*Buoy_R^3*2*zetaot.*zetaotx;

F_viscous = 0.5*rho*2*Buoy_R*Cd*trapz(u1.*abs(u1),2)'*dz...
                                    +0.5*rho*2*Buoy_R*Cd*zetao.*(u1(:,1).*abs(u1(:,1)))';
%% Plots
FIT = F_Morison+F_body_nonlinear+F_wave_ele_nonl+F_crest+F_second;   %without viscous force 
% Compare FIT and experiment 
figure;
plot(time_record,FIT,'b');
hold on;
plot(time_record,Fx,'r');
%title('FIT (small ka) vs Experiment');
legend('FIT','Experiment');
xlabel('Time [s]');
ylabel('Loads [N]');
set(gca,'FontSize',20,'FontWeight','bold');
% Compare FIT, finite Ka and experiment (first 100 seconds)
figure;
plot(time_record(1:1000),-(F_FK_SB(1,:)+F_FK_DS(1,:)+F_FK_SW(1,:)+F_D_SB(1,:)+F_D_DS(1,:)+F_H_SB(1,:)+F_H_SW(1,:)+F_FS_SB'),'b');
hold on;
plot(time_record(1:1000),FIT(1:1000),'g');
hold on;
plot(time_record(1:1000),F_second(1:1000),'y');
hold on;
plot(time_record(1:1000),Fx(1:1000),'r');
%title('FIT (small ka) vs Experiment');
legend('FIT (finte Ka formula with free surface force)','FIT(small Ka with second-order force)','Second-order force','Experiment');
xlabel('Time [s]');
ylabel('Loads [N]');
set(gca,'FontSize',20,'FontWeight','bold');

% Plot zeta, zeta_t, zeta_x, zeta_xt with different frequency limit 
figure;
plot(time_record(1:1000),wave(1:1000));
hold on;
plot(time_record(1:1000),zetao(1:1000));
legend('Experimental wave elevation','Simulated wave elevation with frequency limit 2 rad/s');
xlabel('Time [sec]');
ylabel('Wave [m]');
set(gca,'FontSize',20,'FontWeight','bold');
% Compare FIT and experiment 
figure;
plot(time_record(1:1000),F_Morison(1:1000),'b');
hold on;
plot(time_record(1:1000),FIT(1:1000));
hold on;
plot(time_record(1:1000),Fx(1:1000),'r');
%title('FIT (small ka) vs Experiment');
legend('FIT linear force with frequency limit 2 rad/s','FIT total simulated force','Experimental force');
xlabel('Time [s]');
ylabel('Loads [N]');
set(gca,'FontSize',20,'FontWeight','bold');



% Load F_FS_DD and plot F_FS_DD
load('F_FS_DD_1_2');
figure;
plot(time_record(1:1000),-(F_FK_SB(1,:)+F_FK_DS(1,:)+F_FK_SW(1,:)+F_D_SB(1,:)+F_D_DS(1,:)+F_H_SB(1,:)+F_H_SW(1,:)+F_FS_SB'+F_FS_DD_1_2),'b');
hold on;
plot(time_record(1:1000),FIT(1:1000),'g');
hold on;
plot(time_record(1:1000),Fx(1:1000),'r');
hold on;
plot(time_record(1:1000),F_FS_DD_1_2(1:1000),'color',[0,0,0]);
%title('FIT (small ka) vs Experiment');
legend('FIT with free surface impulse(ID and DD)','FIT (small Ka formula)','Experiment','F\_FS\_DD');
xlabel('Time [s]');
ylabel('Loads [N]');
set(gca,'FontSize',20,'FontWeight','bold');
% figure;
% plot(time_record,F_viscous,'b');
% hold on;
% plot(time_record,Fx,'r');
% %title('Viscous Force vs. Experiment');
% legend('Viscous','Experiment');
% xlabel('Time [s]');
% ylabel('Loads [N]');
% set(gca,'FontSize',20,'FontWeight','bold');

%% Peaks ratio 
ratio = FIT./Fx';
ratiop = ratio(abs(ratio)<2);
% ratio_mean_raw(i) = mean(ratio);
% ratio_mean(i) = mean(ratiop);

figure;
plot(ratio,'b');
%title('FIT (small ka) vs Experiment');
legend('FIT (small Ka assumption)/Experiment');
xlabel('Time [s]');
ylabel( 'Original Ratio Factor');
set(gca,'FontSize',20,'FontWeight','bold');

figure;
plot(ratiop,'b');
%title('FIT (small ka) vs Experiment');
legend('FIT (small Ka assumption)/Experiment');
xlabel('Time [s]');
ylabel( 'Singularity Removed Ratio Factor');
set(gca,'FontSize',20,'FontWeight','bold');

pk1 = findpeaks(FIT);
pk2 = findpeaks(Fx);
pk2 = sort(pk2);
pk2 = pk2(1:length(pk1));
ratio_peak = pk1./pk2';
%% RMS error
rms_force = sqrt(sum((FIT(1:1000)'-Fx(1:1000)).^2)/length(Fx(1:1000)));
rms_wave = sqrt(sum((zetao(1:1000)'-wave(1:1000)).^2)/length(wave(1:1000)));




