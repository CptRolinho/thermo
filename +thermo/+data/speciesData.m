function par = speciesData(speciesName)
%thermo.data.speciesData  Retrieve thermodynamic species metadata (+ optional payload)
%
%   par = thermo.data.speciesData(speciesName)
%
% Inputs
%   speciesName : char/string scalar (case-insensitive match against SpeciesData.Species)
%
% Outputs (struct)
%   par.model, par.formula, par.substance, par.info
%   par.Mm, par.R, par.T_c, par.P_c, par.V_c
%   par.species : canonical key as stored in database
%   par.payload : optional species-specific data block if present in thermo.mat
%
% Notes
% - Robust to working directory changes (uses absolute path to thermo.mat).
% - Uses a persistent cache to avoid repeated disk I/O.
%
% Units expected in database:
%   T [K], P [kPa], V [m^3/kmol], Mm [kg/kmol], R [kJ/kg/K]

    arguments
        speciesName (1,1) {mustBeTextScalar}
    end
    speciesName = string(speciesName);

    % --- Load database (cached) ---
    [SpeciesData, payloadMap] = iLoadDatabase();

    % --- Find species row (case-insensitive) ---
    speciesCol = SpeciesData.Species;
    if iscell(speciesCol)
        speciesCol = string(speciesCol);
    end

    idx = strcmpi(speciesCol, speciesName);
    if ~any(idx)
        avail = unique(speciesCol);
        % Show up to 15 options to keep message readable
        nShow = min(numel(avail), 15);
        msg = sprintf('Unknown species "%s". Example available species: %s%s', ...
            speciesName, strjoin(avail(1:nShow), ", "), ...
            ternary(numel(avail) > nShow, ", ...", ""));
        error('thermo:data:UnknownSpecies', '%s', msg);
    end
    if nnz(idx) > 1
        error('thermo:data:AmbiguousSpecies', ...
            'Species "%s" matched %d rows. Database must have unique Species keys.', ...
            speciesName, nnz(idx));
    end

    row = SpeciesData(idx, :);

    % --- Build output struct (normalize types) ---
    par = struct();
    par.species   = string(row.Species{1});      % canonical key
    par.model     = string(row.Model{1});
    par.formula   = string(row.Formula{1});
    par.substance = string(row.Substance{1});
    par.info      = string(row.Information{1});

    par.Mm  = iScalar(row.M);     % kg/kmol
    par.R   = iScalar(row.R);     % kJ/kg/K
    par.T_c = iScalar(row.T_c);   % K
    par.P_c = iScalar(row.P_c);   % kPa
    par.V_c = iScalar(row.V_c);   % m^3/kmol

    % --- Optional payload (e.g., Air, R134a, Water cells/structs used by models) ---
    par.payload = [];
    if ~isempty(payloadMap) && isKey(payloadMap, lower(par.species))
        par.payload = payloadMap(lower(par.species));
    end
end

% ========================= Helpers =========================

function [SpeciesData, payloadMap] = iLoadDatabase()
    persistent cachedSpeciesData cachedPayloadMap cachedMatPath

    if ~isempty(cachedSpeciesData)
        SpeciesData = cachedSpeciesData;
        payloadMap  = cachedPayloadMap;
        return
    end

    % Resolve absolute path to thermo.mat based on this file location
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    matPath  = fullfile(thisDir, 'thermo.mat');

    if exist(matPath, 'file') ~= 2
        error('thermo:data:MissingDatabase', ...
            'Could not find thermo.mat at: %s', matPath);
    end

    S = load(matPath); % load all to discover available payload variables

    if ~isfield(S, 'SpeciesData')
        error('thermo:data:InvalidDatabase', ...
            'thermo.mat must contain variable "SpeciesData" (table). Found: %s', ...
            strjoin(fieldnames(S), ', '));
    end

    SpeciesData = S.SpeciesData;
    if ~istable(SpeciesData)
        error('thermo:data:InvalidDatabase', ...
            '"SpeciesData" must be a table.');
    end

    iValidateSpeciesTable(SpeciesData);

    % Build optional payload map: if thermo.mat contains variables named
    % like 'Air', 'R134', 'Water' (cells/structs), map them by lowercase key.
    payloadMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    fn = fieldnames(S);
    for k = 1:numel(fn)
        name = fn{k};
        if strcmp(name, 'SpeciesData')
            continue
        end
        % Heuristic: accept payload vars that are cell or struct
        if iscell(S.(name)) || isstruct(S.(name))
            payloadMap(lower(name)) = S.(name);
        end
    end

    % Cache
    cachedSpeciesData = SpeciesData;
    cachedPayloadMap  = payloadMap;
    cachedMatPath     = matPath; %#ok<NASGU>
end

function iValidateSpeciesTable(T)
    required = ["Species","Model","Formula","Substance","Information","M","R","T_c","P_c","V_c"];
    missing = required(~ismember(required, string(T.Properties.VariableNames)));
    if ~isempty(missing)
        error('thermo:data:InvalidDatabase', ...
            'SpeciesData is missing required columns: %s', strjoin(missing, ', '));
    end

    % Basic type sanity checks
    if height(T) < 1
        error('thermo:data:InvalidDatabase', 'SpeciesData table is empty.');
    end

    % Ensure Species keys are nonempty and unique (case-insensitive)
    sp = T.Species;
    if iscell(sp); sp = string(sp); end
    sp = strip(string(sp));

    if any(strlength(sp) == 0)
        error('thermo:data:InvalidDatabase', 'SpeciesData.Species contains empty keys.');
    end

    spLower = lower(sp);
    if numel(unique(spLower)) ~= numel(spLower)
        error('thermo:data:InvalidDatabase', ...
            'SpeciesData.Species contains duplicate keys (case-insensitive).');
    end
end

function y = iScalar(x)
    % Extract numeric scalar from table content robustly
    if iscell(x)
        x = x{1};
    end
    validateattributes(x, {'numeric'}, {'scalar','real','finite'});
    y = double(x);
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
