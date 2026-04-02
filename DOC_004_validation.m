close all, clear all, clc
%% Load data
meas = load("dataset\BID002_RANDCh_30052024.xlsx");
t = meas(:,1)';
y = meas(:,2)';
u = meas(:,3)';
load("DS_002_RCpar.mat")
par = [R0_c R1_c R2_c C1_c C2_c];
Qn = 3.08;
%% Simulation
soc(1) = 1;
for ii = 2:length(t)
    soc(ii) = soc(ii-1) - u(ii)/(3600*Qn);
end
x(:,1) = [0;0;soc(1)];
phi(1) = pVoc*(x(3,1).^(length(pVoc)-1:-1:0))';
ym(1) = [-1 -1 0]*x - R0_c*u(1) + phi(1);
for ii = 2:length(t)
    [x(:,ii),ym(ii),phi(ii)] = ecm2(x(:,ii-1),u(ii),u(ii-1),par,pVoc);
end
%% Plot figures
figure()
plot(t,u,'k-','Linewidth',2)
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',16)
xlabel("Time (s)",'FontSize',16,'Interpreter','latex')
ylabel("Current (A)",'FontSize',16,'Interpreter','latex')
xlim([0 length(t)-1])

figure()
plot(t,soc,'k-','LineWidth',2)
hold on
plot(t,x(3,:),'r--','Linewidth',2)
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',16)
xlabel("Time (s)",'FontSize',16,'Interpreter','latex')
ylabel("SOC",'FontSize',16,'Interpreter','latex')
legend({"Measured","Model"},'FontSize',14,'Interpreter','latex')
xlim([0 length(t)-1])

figure()
plot(t,y,'k-','LineWidth',2)
hold on
plot(t,ym,'r--','Linewidth',2)
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',16)
xlabel("Time (s)",'FontSize',16,'Interpreter','latex')
ylabel("Voltage (V)",'FontSize',16,'Interpreter','latex')
legend({"Measured","Model"},'FontSize',14,'Interpreter','latex')
xlim([0 length(t)-1])
ylim([2.5 4.2])