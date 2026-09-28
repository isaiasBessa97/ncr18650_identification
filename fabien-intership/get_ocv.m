function [soc_axe, V_average, Qn] = get_ocv(file_charge, file_discharge)
% GET_PSEUDO_OCV Calculates the average OCV curve from charge and discharge files
%
% Inputs:
%   file_charge    : Full path to the charge file (string)
%   file_discharge : Full path to the discharge file (string)
%
% Outputs:
%   soc_axe        : Common SoC axis (0 to 100%)
%   V_average      : Resulting Pseudo-OCV voltage vector
%   Qn             : Calculated nominal capacity (Ah)

    %% 1. Loading Data (Added robustness)
    % Replacing the old readmatrix calls with a custom robust loading function
    data_ch = load_battery_data(file_charge);
    data_dis = load_battery_data(file_discharge);

    % Extract columns
    V_ch = data_ch(:, 2);  I_ch = abs(data_ch(:, 3));
    V_dis = data_dis(:, 2); I_dis = abs(data_dis(:, 3));

    % Calculate nominal capacity (Qn) in Ah based on discharge
    Qn = sum(I_dis) / 3600;

    %% 2. SoC Calculation (Coulomb Counting)
    % --- CHARGE ---
    soc_ch = zeros(length(V_ch), 1);
    soc_ch(1) = 0;
    for ii = 2:length(V_ch)
        soc_ch(ii) = soc_ch(ii-1) + (1 / (3600 * Qn)) * I_ch(ii);
    end

    % --- DISCHARGE ---
    soc_dis = zeros(length(V_dis), 1);
    soc_dis(1) = 1;
    for ii = 2:length(V_dis)
        soc_dis(ii) = soc_dis(ii-1) - (1 / (3600 * Qn)) * I_dis(ii);
    end

    soc_ch_pct = soc_ch * 100;
    soc_dis_pct = soc_dis * 100;

    %% 3. Interpolation and Average
    soc_axe = linspace(0, 100, 1000)';
    V_ch_aligned = interp1(soc_ch_pct, V_ch, soc_axe, 'linear');
    V_dis_aligned = interp1(soc_dis_pct, V_dis, soc_axe, 'linear');

    V_average = (V_ch_aligned + V_dis_aligned) / 2;
end

%% --- Helper function for robust file reading ---
function data = load_battery_data(filename)
    % Use fileparts to extract the file extension
    [~, ~, ext] = fileparts(filename);
    
    if strcmpi(ext, '.txt')
        % Original strict behavior: ensures full backward compatibility with old scripts
        data = readmatrix(filename, 'Delimiter', ';', 'NumHeaderLines', 1);
        
    elseif strcmpi(ext, '.csv')
        % For CSV files, let MATLAB automatically detect the delimiter (, or ;)
        try
            % Ignore the first line assuming it is a header
            data = readmatrix(filename, 'NumHeaderLines', 1);
        catch
            % Robust fallback option in case the CSV has a complex structure
            opts = detectImportOptions(filename);
            data = readmatrix(filename, opts);
        end
        
    else
        % Fallback for any other file extension
        data = readmatrix(filename);
    end
end