close all, clear all, clc
meas = readtable("dataset-thermal-sem_rele\BID003\BID003_MPDch_17042026.txt");
time = meas.time;
volt = meas.voltage;
curr = meas.current;
temp_s = meas.surface_temperature_plus;
temp_a = meas.ambient_temperature;

figure()
plot(time,curr)
xlabel("Time (s)")
ylabel("Current (A)")

figure()
plot(time,volt)
xlabel("Time (s)")
ylabel("Voltage (V)")
ylim([2.5 4.2])

figure()
plot(time,temp_s)
xlabel("Time (s)")
ylabel("Surface temperature (ºC)")

figure()
plot(time,temp_a)
xlabel("Time (s)")
ylabel("Ambient temperature (ºC)")