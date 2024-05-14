function par = speciesData(speciesName)
%SPECIES Thermodynamic species database
%   Retrive species model, formula, substance, molar mass, gas constant,
%   critical-point temperature, pressure and volume.
% Load species data table
    load('thermo.mat', 'SpeciesData');

    % Find the row in the table
    idx = strcmpi(SpeciesData.Species, speciesName);
    if any(idx)
        % Extract data for the species
        rowData = SpeciesData(idx, :);        
        par.model = rowData.Model{1};
        par.formula = rowData.Formula{1};
        par.substance = rowData.Substance{1};
        par.Mm = rowData.M;
        par.R = rowData.R;
        par.T_c = rowData.T_c;
        par.P_c = rowData.P_c;
        par.V_c = rowData.V_c;
    else
        error('Unknown species: %s', speciesName);
    end
end