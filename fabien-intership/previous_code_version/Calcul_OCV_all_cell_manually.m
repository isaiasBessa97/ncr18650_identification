close all; clear all; clc;

% Define file paths 
cell1_charge = 'C:\ncr18650_identification\dataset-thermal\BID001\BID001_CCCV005.0_160102024.txt';
cell1_discharge = 'C:\ncr18650_identification\dataset-thermal\BID001\BID001_CDch005.0_150102024.txt';

cell2_charge = 'C:\ncr18650_identification\dataset-thermal\BID002\BID002_CCCV005.0_19022025.txt';
cell2_discharge = 'C:\ncr18650_identification\dataset-thermal\BID002\BID002_CDch005.0_19022025.txt';

% Call the function for Cell 1
[soc1, ocv1, qn1] = get_ocv(cell1_charge, cell1_discharge);

% Call the function for Cell 2
[soc2, ocv2, qn2] = get_ocv(cell2_charge, cell2_discharge);

% compare them on a single plot
figure; hold on; grid on;
plot(soc1, ocv1, 'LineWidth', 2, 'DisplayName', ['Cell 1 (Qn=' num2str(qn1) 'Ah)']);
plot(soc2, ocv2, 'LineWidth', 2, 'DisplayName', ['Cell 2 (Qn=' num2str(qn2) 'Ah)']);

xlabel('SoC (%)'); ylabel('OCV (V)');
title('Comparison of OCV Curves for Multiple Cells');
legend;