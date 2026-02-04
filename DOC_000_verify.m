close all, clear all, clc
meas = readtable("dataset-thermal\BID002\BID002_CDch005.0_19022025.txt");
time = meas.time;
volt = meas.voltage;
curr = meas.current;

figure()
plot(time,curr)
xlabel("Time (s)")
ylabel("Current (A)")

figure()
plot(time,volt)
xlabel("Time (s)")
ylabel("Voltage (V)")