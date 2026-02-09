function T = listSpecies()
%thermo.data.listSpecies  Return SpeciesData table from thermo.mat (robust path)

    thisFile = mfilename('fullpath');
    dataDir = fileparts(thisFile);
    matPath = fullfile(dataDir, 'thermo.mat');

    if exist(matPath,'file') ~= 2
        error('thermo:data:MissingDatabase', 'Could not find thermo.mat at: %s', matPath);
    end

    S = load(matPath, 'SpeciesData');
    if ~isfield(S,'SpeciesData') || ~istable(S.SpeciesData)
        error('thermo:data:InvalidDatabase', 'thermo.mat must contain a table variable "SpeciesData".');
    end
    T = S.SpeciesData;
end
