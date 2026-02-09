%% smokeTest_Thermo.m
% Smoke test for Thermo class + +thermo package
% Intended for quick validation in student projects and CI-like workflows.
%
% Run:
%   tests/smokeTest_Thermo
%
% Expected:
%   Prints PASS messages; errors on first failure.

clear; clc;

fprintf('=== Thermo smoke test ===\n');

%% 0) Add project root to path (robustly)
thisFile = mfilename('fullpath');
if isempty(thisFile)
    error('SmokeTest:Path', 'Unable to resolve script full path. Save the script and re-run.');
end
testsDir = fileparts(thisFile);
projectRoot = fileparts(testsDir);

addpath(genpath(projectRoot));
fprintf('[INFO] Added to path: %s\n', projectRoot);

%% 1) Resolve key entry points
mustExist('Thermo',       'class');
mustExist('thermo.data.speciesData', 'function');
mustExist('thermo.models.liquid',    'function');
mustExist('thermo.models.idealGas',  'function');

fprintf('[PASS] Entry points resolved.\n');

%% 2) Database / species listing
T = thermo.data.listSpecies();
assert(istable(T) && height(T) >= 1, 'SmokeTest:SpeciesTable', 'SpeciesData table not found or empty.');
req = ["Air","Water","R134a"];
for s = req
    assert(any(strcmpi(string(T.Species), s)), 'SmokeTest:SpeciesMissing', 'Species missing from database: %s', s);
end
fprintf('[PASS] Database loaded; required species present.\n');

%% 3) Ideal gas (class + package)
T0 = 300;           % K
P0 = 101.325;       % kPa
air = Thermo("Air", 'T', T0, 'P', P0);

assertFinite(air.v, 'Air.v');
assertFinite(air.cp, 'Air.cp');
assertFinite(air.cv, 'Air.cv');
assertFinite(air.k, 'Air.k');

vCalc = air.R * air.T / air.P; % (kJ/kg/K * K)/kPa = (kPa*m^3/kg)/kPa = m^3/kg
assertRelClose(air.v, vCalc, 1e-9, 'Ideal gas v consistency');

% Package API
ig = thermo.models.idealGas("Air", T0);
assertFinite(ig.cp, 'idealGas.cp');
fprintf('[PASS] Ideal gas (Air) validated (class + package).\n');

%% 4) Water saturated mixture at 100C
Tw = 373.15;   % K
xw = 0.5;
thW = Thermo("Water", 'T', Tw, 'x', xw);

assert(strcmpi(thW.state, "saturated mixture"), 'SmokeTest:WaterState', 'Expected saturated mixture for Water at T=373.15 K, x=0.5.');

% Sanity: saturation pressure near 1 atm (table-dependent)
assert(abs(thW.P - 101.325) < 2.0, 'SmokeTest:WaterPsat', ...
    'Water saturation pressure at 373.15 K seems off: P=%.6g kPa', thW.P);

% Mixture identities
assertRelClose(thW.u, thW.u_f + xw*thW.u_fg, 1e-9, 'Water mixture u');
assertRelClose(thW.h, thW.h_f + xw*thW.h_fg, 1e-9, 'Water mixture h');
assertRelClose(thW.s, thW.s_f + xw*thW.s_fg, 1e-9, 'Water mixture s');
assertRelClose(thW.v, (1-xw)*thW.v_f + xw*thW.v_g, 1e-9, 'Water mixture v');

fprintf('[PASS] Water saturated mixture validated.\n');

%% 5) Water superheated check (should not be saturated mixture)
% Choose a point above Tsat at 200 kPa (Tsat ~ 393K; choose 500K)
thW2 = Thermo("Water", 'P', 200, 'T', 500);
assert(strlength(thW2.state) > 0, 'SmokeTest:WaterState2', 'Water state should not be empty.');
fprintf('[PASS] Water off-saturation state resolved: %s\n', thW2.state);

%% 6) R134a saturated mixture sanity (P,x)
Pr = 200;   % kPa (must be within your table)
xr = 0.2;
thR = Thermo("R134a", 'P', Pr, 'x', xr);

assert(strcmpi(thR.state, "saturated mixture"), 'SmokeTest:R134aState', ...
    'Expected saturated mixture for R134a at P=200 kPa, x=0.2.');

% Mixture identities (if saturation fields are populated)
assertFinite(thR.u_f, 'R134a.u_f');
assertFinite(thR.u_fg,'R134a.u_fg');
assertRelClose(thR.u, thR.u_f + xr*thR.u_fg, 1e-9, 'R134a mixture u');
assertRelClose(thR.h, thR.h_f + xr*thR.h_fg, 1e-9, 'R134a mixture h');
assertRelClose(thR.s, thR.s_f + xr*thR.s_fg, 1e-9, 'R134a mixture s');
assertRelClose(thR.v, (1-xr)*thR.v_f + xr*thR.v_g, 1e-9, 'R134a mixture v');

fprintf('[PASS] R134a saturated mixture validated.\n');

%% 7) Package call equivalence (spot-check)
pw = thermo.models.liquid("Water", 'T', Tw, 'x', xw);
assertRelClose(pw.P, thW.P, 1e-12, 'Class vs package P (Water T,x)');
assertRelClose(pw.h, thW.h, 1e-12, 'Class vs package h (Water T,x)');
fprintf('[PASS] Class vs package equivalence spot-check.\n');

fprintf('\n=== ALL SMOKE TESTS PASSED ✅ ===\n');

%% ---------- Local helper functions ----------
function mustExist(name, kind)
    w = which(name);
    if isempty(w)
        error('SmokeTest:Missing', 'Missing %s: %s', kind, name);
    end
end

function assertFinite(x, label)
    if ~(isscalar(x) && isnumeric(x) && isfinite(x))
        error('SmokeTest:NonFinite', '%s must be finite scalar numeric. Got: %s', label, mat2str(x));
    end
end

function assertRelClose(a, b, tol, label)
    if ~(isscalar(a) && isscalar(b) && isnumeric(a) && isnumeric(b))
        error('SmokeTest:BadScalar', '%s inputs must be scalar numerics.', label);
    end
    denom = max(1, abs(b));
    rel = abs(a - b) / denom;
    if rel > tol
        error('SmokeTest:NotClose', '%s failed: a=%.12g, b=%.12g, rel=%.3g > tol=%.3g', ...
            label, a, b, rel, tol);
    end
end
