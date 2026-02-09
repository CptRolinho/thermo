function setupThermo()
%setupThermo  One-command setup for Thermo (class + package)
%
% Adds the project root to the MATLAB path, checks availability, and runs
% the smoke test if present.
%
% Run:
%   setupThermo

    fprintf('=== Thermo setup ===\n');

    % Resolve this file location -> project root
    thisFile = mfilename('fullpath');
    if isempty(thisFile)
        error('setupThermo:Path', 'Unable to resolve setupThermo full path.');
    end
    projectRoot = fileparts(thisFile);

    % Add ONLY project root (never add +thermo directly)
    addpath(projectRoot);

    % Add common non-namespace folders if present
    addIfExists(fullfile(projectRoot, 'tests'));
    addIfExists(fullfile(projectRoot, 'examples'));

    fprintf('[INFO] Added to path: %s\n', projectRoot);

    % Verify entry points
    mustExist('Thermo', 'class');
    mustExist('thermo.data.speciesData', 'function');
    mustExist('thermo.models.liquid', 'function');
    mustExist('thermo.models.idealGas', 'function');
    fprintf('[PASS] Entry points resolved.\n');

    % Run smoke test if available
    smokeCandidates = {
        fullfile(projectRoot, 'tests', 'smokeTest_Thermo.m')
        fullfile(projectRoot, 'smokeTest_Thermo.m')
        fullfile(projectRoot, '+thermo', 'tests', 'smokeTest_Thermo.m')
        };

    ranSmoke = false;
    for i = 1:numel(smokeCandidates)
        f = smokeCandidates{i};
        if exist(f, 'file') == 2
            fprintf('[INFO] Running smoke test: %s\n', f);
            run(f);
            ranSmoke = true;
            break
        end
    end

    if ~ranSmoke
        fprintf('[WARN] Smoke test not found. (Recommended: tests/smokeTest_Thermo.m)\n');
    end

    fprintf('\n=== Thermo ready ✅ ===\n');
end

function addIfExists(p)
    if exist(p, 'dir') == 7
        addpath(p);
    end
end

function mustExist(name, kind)
    w = which(name);
    if isempty(w)
        error('setupThermo:Missing', 'Missing %s: %s', kind, name);
    end
end
