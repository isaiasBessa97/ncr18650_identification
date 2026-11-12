clear all, close all, clc
%% load
dataMeas = readtable("BID004_NCM_1_MPDCH_20261107_135608.csv");

volt = dataMeas.V;
curr = dataMeas.I;
temp = dataMeas.T_C;
time = 0:1:length(volt)-1;

%% Plot

figure()
plot(time,curr,'LineWidth',2)
grid on
set(gca,'ticklabelinterpreter','latex','fontsize',22)
xlabel("Time (s)","Interpreter","latex",'FontSize',22)
ylabel("Current (A)","Interpreter","latex",'FontSize',22)

figure()
plot(time,volt,'LineWidth',2)
grid on
set(gca,'ticklabelinterpreter','latex','fontsize',22)
xlabel("Time (s)","Interpreter","latex",'FontSize',22)
ylabel("Voltage (V)","Interpreter","latex",'FontSize',22)

figure()
plot(time,temp,'LineWidth',2)
grid on
set(gca,'ticklabelinterpreter','latex','fontsize',22)
xlabel("Time (s)","Interpreter","latex",'FontSize',22)
ylabel("Temperature (C)","Interpreter","latex",'FontSize',22)