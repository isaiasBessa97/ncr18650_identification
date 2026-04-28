function soc_array = get_soc(file_path, Qn, initial_soc)
% GET_SOC Calculates the continuous State of Charge (SoC) for a given test
%
% Inputs:
%   file_path   : Full path to the test file
%   Qn          : Nominal capacity of the cell in Ah
%   initial_soc : Starting SoC of the test (e.g., 100 for a fully charged cell)
%
% Outputs:
%   soc_array   : Array containing the SoC value for every timestamp

    % Read Data
    data = readmatrix(file_path, 'Delimiter', ';', 'NumHeaderLines', 1);
    time = data(:, 1);
    I = data(:, 3); 
    
    % Initialize SoC array
    num_points = length(I);
    soc_array = zeros(num_points, 1);
    soc_array(1) = initial_soc;
    
    % Coulomb Counting loop
    for ii = 2:num_points
        % Calculate dt (Time step in seconds)
        dt = time(ii) - time(ii-1);
        if dt == 0; dt = 1; end % Fallback safety
        
        % Calculate SoC in percentage (added 100 * multiplier)
        soc_array(ii) = soc_array(ii-1) - (100 * dt / (3600 * Qn)) * I(ii);
    end
end