function par = liquid(species, varargin)
%thermo.models.liquid  Generic liquid-vapor thermo model (Water, R134a, ...)
% Reference tables: Çengel (as stored in thermo.mat payloads)
%
% Usage:
%   par = thermo.models.liquid("Water",'T',373.15,'x',0.5)
%   par = thermo.models.liquid("R134a",'P',200,'T',310)
%
% Inputs (Name-Value):
%   T [K], P [kPa], x [-], u [kJ/kg], v [m^3/kg], h [kJ/kg], s [kJ/kg/K]
%
% Output:
%   par struct with fields: T,P,v,u,h,s,x,state and saturation fields.

    arguments
    species (1,1) {mustBeTextScalar}
    end
    arguments (Repeating)
        varargin
    end


    % ---- parse inputs ----
    p = inputParser;
    p.FunctionName = 'thermo.models.liquid';

    % validation (allow NaN to mean "not provided")
    vPos = @(x) (isnan(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
    vAny = @(x) (isnan(x) || (isnumeric(x) && isscalar(x) && isfinite(x)));
    vX   = @(x) (isnan(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x>=0 && x<=1));

    addParameter(p,'T',NaN,vPos);
    addParameter(p,'P',NaN,vPos);
    addParameter(p,'x',NaN,vX);
    addParameter(p,'u',NaN,vAny);
    addParameter(p,'v',NaN,vPos);
    addParameter(p,'h',NaN,vAny);
    addParameter(p,'s',NaN,vAny);

    parse(p,varargin{:});

    par = initPar(p.Results);

    provided = ~isnan([par.T, par.P, par.x, par.u, par.v, par.h, par.s]);
    if sum(provided) < 2
        error('thermo:models:InvalidState', ...
            'Provide at least two among {T,P,x,u,v,h,s}.');
    end
    if ~isnan(par.x) && isnan(par.T) && isnan(par.P)
        error('thermo:models:InvalidState', ...
            'If quality x is provided, you must also provide T or P.');
    end

    % ---- load payload (cached) ----
    payload = thermo.data.getPayload(species);

    tblT = payload{1};
    tblP = payload{2};

    % tolerances
    tolT = @(Tsat) 1e-6 * max(1, abs(Tsat));
    tolP = @(P)    1e-8 * max(1, abs(P)); % pressure grids usually exact

    % ---- main logic ----
    if ~isnan(par.T) && ~isnan(par.P)
        Tsat = interpP2Tsat(tblP, par.P);
        par.T_sat = Tsat;

        if par.T < Tsat - tolT(Tsat)
            par.state = "compressed liquid";
            par = compressedRegion(par, payload, tolP);
        elseif abs(par.T - Tsat) <= tolT(Tsat)
            par.state = "saturated mixture";
            par = saturatedFromP(par, tblP);
            par.T = par.T_sat;
            par = resolveMixture(par);
        else
            par.state = "superheated vapor";
            par = superheatedRegion(par, payload, tolP);
        end

    elseif ~isnan(par.T) && ~isnan(par.x)
        par.state = "saturated mixture";
        par = saturatedFromT(par, tblT);
        par.P = par.P_sat;
        par = resolveMixture(par);

    elseif ~isnan(par.P) && ~isnan(par.x)
        par.state = "saturated mixture";
        par = saturatedFromP(par, tblP);
        par.T = par.T_sat;
        par = resolveMixture(par);

    else
        yProp = firstProvidedUVHS(par);

        if isnan(yProp)
            error('thermo:models:InvalidState', ...
                'Unsupported combination. Provide (T,P) or include x.');
        end

        if isnan(par.T) && isnan(par.P)
            error('thermo:models:InvalidState', ...
                'When providing u/v/h/s, also provide T or P.');
        end

        if ~isnan(par.T)
            sat = saturationThresholdsFromT(tblT, par.T, yProp);
            y = par.(yProp);

            if y < sat.yf
                par.state = "compressed liquid";
                if isnan(par.P), par.P = sat.Psat; end
                par = compressedRegion(par, payload, tolP);
            elseif y < sat.yg
                par.state = "saturated mixture";
                par = saturatedFromT(par, tblT);
                par.P = par.P_sat;
                par.x = (y - sat.yf) / sat.yfg;
                par = resolveMixture(par);
            else
                par.state = "superheated vapor";
                if isnan(par.P), par.P = sat.Psat; end
                par = superheatedRegion(par, payload, tolP);
            end

        else
            sat = saturationThresholdsFromP(tblP, par.P, yProp);
            y = par.(yProp);

            if y < sat.yf
                par.state = "compressed liquid";
                par.T = sat.Tsat;
                par = compressedRegion(par, payload, tolP);
            elseif y < sat.yg
                par.state = "saturated mixture";
                par = saturatedFromP(par, tblP);
                par.T = par.T_sat;
                par.x = (y - sat.yf) / sat.yfg;
                par = resolveMixture(par);
            else
                par.state = "superheated vapor";
                par.T = sat.Tsat;
                par = superheatedRegion(par, payload, tolP);
            end
        end
    end
end

% ========================= Local helpers =========================

function par = initPar(r)
    par = struct();
    par.T = r.T; par.P = r.P; par.x = r.x;
    par.u = r.u; par.v = r.v; par.h = r.h; par.s = r.s;

    par.T_sat=[]; par.P_sat=[];
    par.v_f=[]; par.v_g=[];
    par.u_f=[]; par.u_g=[]; par.u_fg=[];
    par.h_f=[]; par.h_g=[]; par.h_fg=[];
    par.s_f=[]; par.s_g=[]; par.s_fg=[];
    par.u_i=[]; par.u_ig=[];
    par.h_i=[]; par.h_ig=[];
    par.s_i=[]; par.s_ig=[];
    par.state="";
end

function Tsat = interpP2Tsat(tblP, P)
    Tsat = interp1safe(tblP.P, tblP.Tsat, P);
end

function Psat = interpT2Psat(tblT, T)
    Psat = interp1safe(tblT.T, tblT.Psat, T);
end

function par = saturatedFromT(par, tblT)
    par.T_sat = par.T;
    par.P_sat = interpT2Psat(tblT, par.T);
    par.v_f   = interp1safe(tblT.T, tblT.vf,  par.T);
    par.v_g   = interp1safe(tblT.T, tblT.vg,  par.T);
    par.u_f   = interp1safe(tblT.T, tblT.uf,  par.T);
    par.u_fg  = interp1safe(tblT.T, tblT.ufg, par.T);
    par.u_g   = interp1safe(tblT.T, tblT.ug,  par.T);
    par.h_f   = interp1safe(tblT.T, tblT.hf,  par.T);
    par.h_fg  = interp1safe(tblT.T, tblT.hfg, par.T);
    par.h_g   = interp1safe(tblT.T, tblT.hg,  par.T);
    par.s_f   = interp1safe(tblT.T, tblT.sf,  par.T);
    par.s_fg  = interp1safe(tblT.T, tblT.sfg, par.T);
    par.s_g   = interp1safe(tblT.T, tblT.sg,  par.T);
end

function par = saturatedFromP(par, tblP)
    par.P_sat = par.P;
    par.T_sat = interp1safe(tblP.P, tblP.Tsat, par.P);
    par.v_f   = interp1safe(tblP.P, tblP.vf,   par.P);
    par.v_g   = interp1safe(tblP.P, tblP.vg,   par.P);
    par.u_f   = interp1safe(tblP.P, tblP.uf,   par.P);
    par.u_fg  = interp1safe(tblP.P, tblP.ufg,  par.P);
    par.u_g   = interp1safe(tblP.P, tblP.ug,   par.P);
    par.h_f   = interp1safe(tblP.P, tblP.hf,   par.P);
    par.h_fg  = interp1safe(tblP.P, tblP.hfg,  par.P);
    par.h_g   = interp1safe(tblP.P, tblP.hg,   par.P);
    par.s_f   = interp1safe(tblP.P, tblP.sf,   par.P);
    par.s_fg  = interp1safe(tblP.P, tblP.sfg,  par.P);
    par.s_g   = interp1safe(tblP.P, tblP.sg,   par.P);
end

function par = resolveMixture(par)
    if isnan(par.x)
        yProp = firstProvidedUVHS(par);
        if isnan(yProp)
            warning('thermo:models:MixtureUnderdetermined', ...
                'Saturated mixture detected but insufficient data to compute quality x.');
            return
        end
        y = par.(yProp);
        yf  = par.(strcat(yProp,'_f'));
        yfg = par.(strcat(yProp,'_fg'));
        par.x = (y - yf) / yfg;
    end

    par.x = min(max(par.x, 0), 1);

    par.u = par.u_f + par.x * par.u_fg;
    par.h = par.h_f + par.x * par.h_fg;
    par.s = par.s_f + par.x * par.s_fg;
    par.v = par.x * par.v_g + (1 - par.x) * par.v_f;
end

function yProp = firstProvidedUVHS(par)
    if ~isnan(par.u); yProp = "u"; return; end
    if ~isnan(par.v); yProp = "v"; return; end
    if ~isnan(par.h); yProp = "h"; return; end
    if ~isnan(par.s); yProp = "s"; return; end
    yProp = NaN;
end

function sat = saturationThresholdsFromT(tblT, T, yProp)
    sat.Psat = interpT2Psat(tblT, T);
    sat.yf  = interp1safe(tblT.T, tblT.(yProp+"f"),  T);
    sat.yg  = interp1safe(tblT.T, tblT.(yProp+"g"),  T);
    sat.yfg = interp1safe(tblT.T, tblT.(yProp+"fg"), T);
end

function sat = saturationThresholdsFromP(tblP, P, yProp)
    sat.Tsat = interp1safe(tblP.P, tblP.Tsat, P);
    sat.yf  = interp1safe(tblP.P, tblP.(yProp+"f"),  P);
    sat.yg  = interp1safe(tblP.P, tblP.(yProp+"g"),  P);
    sat.yfg = interp1safe(tblP.P, tblP.(yProp+"fg"), P);
end

function par = superheatedRegion(par, payload, tolP)
    if numel(payload) < 3 || isempty(payload{3})
        error('thermo:models:MissingSuperheatedTables', ...
            'Payload has no superheated vapor tables (payload{3}).');
    end
    tabCell = payload{3};
    Pgrid = cell2mat(tabCell(:,1));
    idx = findNearestIdx(Pgrid, par.P, tolP(par.P));

    par.T_sat = tabCell{idx,2};
    tab = tabCell{idx,3};

    par.v = interp1safe(tab.T, tab.v, par.T);
    par.u = interp1safe(tab.T, tab.u, par.T);
    par.h = interp1safe(tab.T, tab.h, par.T);
    par.s = interp1safe(tab.T, tab.s, par.T);

    par.x = [];
end

function par = compressedRegion(par, payload, tolP)
    tblT = payload{1};

    if numel(payload) >= 4 && ~isempty(payload{4})
        compCell = payload{4};
        Pmin = compCell{1,1};

        if par.P < Pmin
            par = saturatedFromT(par, tblT);
            par.v_g=[]; par.u_fg=[]; par.u_g=[]; par.h_fg=[]; par.h_g=[]; par.s_fg=[]; par.s_g=[];
            par.v = par.v_f; par.u = par.u_f; par.h = par.h_f; par.s = par.s_f;
            par.x = [];
            return
        end

        Pgrid = cell2mat(compCell(:,1));
        idx = findNearestIdx(Pgrid, par.P, tolP(par.P));

        par.T_sat = compCell{idx,2};
        tab = compCell{idx,3};

        par.v = interp1safe(tab.T, tab.v, par.T);
        par.u = interp1safe(tab.T, tab.u, par.T);
        par.h = interp1safe(tab.T, tab.h, par.T);
        par.s = interp1safe(tab.T, tab.s, par.T);

        par.x = [];
    else
        par = saturatedFromT(par, tblT);
        par.v_g=[]; par.u_fg=[]; par.u_g=[]; par.h_fg=[]; par.h_g=[]; par.s_fg=[]; par.s_g=[];
        par.v = par.v_f; par.u = par.u_f; par.h = par.h_f; par.s = par.s_f;
        par.x = [];
    end
end

function idx = findNearestIdx(grid, q, tol)
    [dmin, idx] = min(abs(grid - q));
    if dmin > tol
        error('thermo:models:PressureNotFound', ...
            'Pressure %.8g kPa not found in grid (min diff %.3g > tol %.3g).', q, dmin, tol);
    end
end
