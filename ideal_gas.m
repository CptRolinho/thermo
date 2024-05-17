function par = ideal_gas(species, T)
    % Parameters for thermodynamic model for a pseudo species as an ideal gas
    % Reference: Çengel
    % v.0.2.0

    % The data file which contains all species data
    fileName = 'thermo.mat';

    % Dynamically access the species data from the intended file
    % Perform a safety verification and access
    speciesData = safelyLoadVariable(fileName, species);

    % Interpolation to find properties at the given temperature T
    % Using first cell element of speciesData for cv, cp, and k,
    % and second cell element for h, u, and s_0.
    par.cv = interp1(speciesData{1}.T, speciesData{1}.cv, T);   % Specific heat capacity at constant volume (kJ/kg*K)
    par.cp = interp1(speciesData{1}.T, speciesData{1}.cp, T);   % Specific heat capacity at constant pressure (kJ/kg*K)
    par.k = interp1(speciesData{1}.T, speciesData{1}.k, T);     % Adiabatic index
    par.h = interp1(speciesData{2}.T, speciesData{2}.h, T);     % Enthalpy (kJ/kg)
    par.u = interp1(speciesData{2}.T, speciesData{2}.u, T);     % Internal energy (kJ/kg)
    par.s_0 = interp1(speciesData{2}.T, speciesData{2}.s_0, T); % Entropy function zero (kJ/kg*K)
end

function exists = checkVariableInFile(fileName, varName)
    % Check if a variable exists in a .mat file
    vars = whos('-file', fileName);
    exists = any(strcmp({vars.name}, varName));
end

function data = safelyLoadVariable(fileName, varName)
    % Load a specific variable from a .mat file if it exists
    if checkVariableInFile(fileName, varName)
        S = load(fileName, varName);
        data = S.(varName);
    else
        error('Data on species not found.');
    end
end
