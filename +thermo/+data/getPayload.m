function payload = getPayload(species)
%thermo.data.getPayload  Return species payload stored in thermo.mat (cached)
%
% payload = thermo.data.getPayload("Water")
% payload is whatever is stored as variable Water/Air/R134a in thermo.mat (cell/struct)

    arguments
        species (1,1) {mustBeTextScalar}
    end
    species = string(species);

    persistent cacheMap
    if isempty(cacheMap)
        cacheMap = iLoadAllPayloads();
    end

    key = lower(species);
    if ~isKey(cacheMap, key)
        error('thermo:data:MissingPayload', ...
            'No payload variable "%s" found in thermo.mat.', species);
    end
    payload = cacheMap(key);
end

function map = iLoadAllPayloads()
    matPath = iMatPath();
    S = load(matPath); % load all; small file, and cached once
    fn = fieldnames(S);

    map = containers.Map('KeyType','char','ValueType','any');
    for k = 1:numel(fn)
        name = fn{k};
        if strcmp(name,'SpeciesData')
            continue
        end
        if iscell(S.(name)) || isstruct(S.(name))
            map(lower(name)) = S.(name);
        end
    end
end

function matPath = iMatPath()
    % thermo.mat is located next to this file
    thisFile = mfilename('fullpath');
    thisDir = fileparts(thisFile);
    matPath = fullfile(thisDir, 'thermo.mat');
    if exist(matPath,'file') ~= 2
        error('thermo:data:MissingDatabase', 'Could not find thermo.mat at: %s', matPath);
    end
end
