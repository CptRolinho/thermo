function par = liquid_Water(varargin)
%Parameters for thermodynamic model for air as a pseudo species
%Ref: Çengel
load('thermo.mat', 'Water');

p = inputParser;
addParameter(p, 'T', NaN);
addParameter(p, 'P', NaN);
addParameter(p, 'x', NaN);
addParameter(p, 'u', NaN);
parse(p, varargin{:});

par.T = p.Results.T;    % Temperature (K)
par.P = p.Results.P;    % Pressure (kPa)
par.x = p.Results.x;    % Quality
par.u = p.Results.u;    % Internal energy (kJ/kg)

par.v     = []; % Specific volume (m3/kg)
par.h     = []; % Enthalpy (kJ/kg)
par.s     = []; % Entropy (kJ/kg*K)
par.T_sat = []; % Saturated temperature (K)
par.P_sat = []; % Saturated pressure (kPa)
par.v_f   = []; % Saturated liquid specific volume (m3/kg)
par.v_g   = []; % Saturated vapor specific volume (m3/kg)
par.u_f   = []; % Saturated liquid internal energy (kJ/kg)
par.u_g   = []; % Saturated vapor internal energy (kJ/kg)
par.u_fg  = []; % Evaporation internal energy (kJ/kg)
par.u_i   = []; % Saturated ice internal energy (kJ/kg)
par.u_ig  = []; % Sublimation internal energy (kJ/kg)
par.h_f   = []; % Saturated liquid enthalpy (kJ/kg)
par.h_g   = []; % Saturated vapor enthalpy (kJ/kg)
par.h_fg  = []; % Evaporation enthalpy (kJ/kg)
par.h_i   = []; % Saturated ice enthalpy (kJ/kg)
par.h_ig  = []; % Sublimation enthalpy (kJ/kg)
par.s_f   = []; % Saturated liquid entropy (kJ/kg*K)
par.s_g   = []; % Saturated vapor entropy (kJ/kg*K)
par.s_fg  = []; % Evaporation entropy (kJ/kg*K)
par.s_i   = []; % Saturated ice entropy (kJ/kg*K)
par.s_ig  = []; % Sublimation entropy (kJ/kg*K)
par.state = []; % Liquid-vapor mixture state

provided = ~isnan([par.T, par.P, par.x, par.u]);
if provided(1) && provided(2) % T and P provided
    par.T_sat = interp1(Water{2}.P,Water{2}.Tsat,par.P);
    if par.T<par.T_sat
        par.state = 'compressed liquid';

        % Get data from compressed liquid
        res=compressed(par,Water);
        par=res;

    elseif par.T==par.T_sat
        par.state = 'saturated mixture';

        % Get data from pressure table
        res=pressure(par,Water);
        par=res;

        par.T = par.T_sat;

        if provided(3)            
            par.u = par.u_f+par.x*par.u_fg;
            par.h = par.h_f+par.x*par.h_fg;
            par.s = par.s_f+par.x*par.s_fg;
            par.v = par.x*par.v_g+(1-par.x)*par.v_f;
        elseif provided(4)
            par.x = (par.u-par.u_f)/par.u_fg;
            par.h = par.h_f+par.x*par.h_fg;
            par.s = par.s_f+par.x*par.s_fg;
            par.v = par.x*par.v_g+(1-par.x)*par.v_f;
        else
            warning('Insufficient data. Specific properties of mixture not calculated')
        end

    else
        par.state = 'superheated vapor';
        % Get data from superheated vapor
        res=superheated(par,Water);
        par=res;
    end

elseif provided(1) && provided(3) % T and x provided
    par.state = 'saturated mixture';

    % Get data from temperature table
    res=temperature(par,Water);
    par=res;

    par.P    = par.P_sat;
    par.u    = par.u_f+par.x*par.u_fg;
    par.h    = par.h_f+par.x*par.h_fg;
    par.s    = par.s_f+par.x*par.s_fg;
    par.v    = par.x*par.v_g+(1-par.x)*par.v_f;

elseif provided(1) && provided(4) % T and u provided
    par.u_f   = interp1(Water{1}.T,Water{1}.uf,par.T);
    par.u_g   = interp1(Water{1}.T,Water{1}.ug,par.T);
    if par.u<par.u_f
        par.state = 'compressed liquid';
        % Get data from compressed liquid
        res=compressed(par,Water);
        par=res;
    elseif par.u<par.u_g
        par.state = 'saturated mixture';

        % Get data from temperature table
        res=temperature(par,Water);
        par=res;

        par.x = (par.u-par.u_f)/par.u_fg;
        par.h = par.h_f+par.x*par.h_fg;
        par.s = par.s_f+par.x*par.s_fg;
        par.v = par.x*par.v_g+(1-par.x)*par.v_f;
    else
        par.state = 'superheated vapor';
        % Get data from superheated vapor
        res=superheated(par,Water);
        par=res;
    end
elseif provided(2) && provided(3) % P and x provided
    par.state = 'saturated mixture';

    % Get data from pressure table
    res=pressure(par,Water);
    par=res;

    par.T    = par.T_sat;
    par.u    = par.u_f+par.x*par.u_fg;
    par.h    = par.h_f+par.x*par.h_fg;
    par.s    = par.s_f+par.x*par.s_fg;
    par.v    = par.x*par.v_g+(1-par.x)*par.v_f;

elseif provided(2) && provided(4) % P and u provided
    par.u_f   = interp1(Water{2}.T,Water{2}.uf,par.P);
    par.u_g   = interp1(Water{2}.T,Water{2}.ug,par.P);

    if par.u<par.u_f
        par.state = 'compressed liquid';
        % Get data from compressed liquid
        res=compressed(par,Water);
        par=res;

    elseif par.u<par.u_g
        par.state = 'saturated mixture';

        % Get data from pressure table
        res=pressure(par,Water);
        par=res;

        par.T = par.T_sat;
        par.x = (par.u-par.u_f)/par.u_fg;
        par.h = par.h_f+par.x*par.h_fg;
        par.s = par.s_f+par.x*par.s_fg;
        par.v = par.x*par.v_g+(1-par.x)*par.v_f;

    else
        par.state = 'superheated vapor';
        % Get data from superheated vapor
        res=superheated(par,Water);
        par=res;
    end
end
end

function res = compressed(par,Water)
res=par;

if par.P < Water{1,4}{1,1}
    % Pressure too small, properties are obtained from Temperature table
    % Get data from temperature table
    res=temperature(par,Water);
    res.v_g   = [];
    res.u_fg  = [];
    res.u_g   = [];
    res.h_fg  = [];
    res.h_g   = [];
    res.s_fg  = [];
    res.s_g   = [];

    res.v = res.v_f;
    res.u = res.u_f;
    res.h = res.h_f;
    res.s = res.s_f;
else
    % Find the row in the table
    idx = find(cell2mat(Water{1,4}(:,1)) == par.P);
    if any(idx)
        % Extract data for the species
        res.T_sat = Water{1,4}{idx,2};

        tab = Water{1,4}{idx,3};
        res.v = interp1(tab.T,tab.v,par.T);
        res.u = interp1(tab.T,tab.u,par.T);
        res.h = interp1(tab.T,tab.h,par.T);
        res.s = interp1(tab.T,tab.s,par.T);
    else
        error('Compressed liquid pressure not found.');
    end
end
res.x=[];
end

function res = pressure(par,Water)
res = par;

res.T_sat = interp1(Water{2}.P,Water{2}.Tsat,par.P);
res.v_f   = interp1(Water{2}.P,Water{2}.vf,par.P);
res.v_g   = interp1(Water{2}.P,Water{2}.vg,par.P);
res.u_f   = interp1(Water{2}.P,Water{2}.uf,par.P);
res.u_fg  = interp1(Water{2}.P,Water{2}.ufg,par.P);
res.u_g   = interp1(Water{2}.P,Water{2}.ug,par.P);
res.h_f   = interp1(Water{2}.P,Water{2}.hf,par.P);
res.h_fg  = interp1(Water{2}.P,Water{2}.hfg,par.P);
res.h_g   = interp1(Water{2}.P,Water{2}.hg,par.P);
res.s_f   = interp1(Water{2}.P,Water{2}.sf,par.P);
res.s_fg  = interp1(Water{2}.P,Water{2}.sfg,par.P);
res.s_g   = interp1(Water{2}.P,Water{2}.sg,par.P);
end

function res = temperature(par,Water)
res = par;

res.P_sat = interp1(Water{1}.T,Water{1}.Psat,par.T);
res.v_f   = interp1(Water{1}.T,Water{1}.vf,par.T);
res.v_g   = interp1(Water{1}.T,Water{1}.vg,par.T);
res.u_f   = interp1(Water{1}.T,Water{1}.uf,par.T);
res.u_fg  = interp1(Water{1}.T,Water{1}.ufg,par.T);
res.u_g   = interp1(Water{1}.T,Water{1}.ug,par.T);
res.h_f   = interp1(Water{1}.T,Water{1}.hf,par.T);
res.h_fg  = interp1(Water{1}.T,Water{1}.hfg,par.T);
res.h_g   = interp1(Water{1}.T,Water{1}.hg,par.T);
res.s_f   = interp1(Water{1}.T,Water{1}.sf,par.T);
res.s_fg  = interp1(Water{1}.T,Water{1}.sfg,par.T);
res.s_g   = interp1(Water{1}.T,Water{1}.sg,par.T);
end

function res=superheated(par,Water)
res = par;

% Find the row in the table
    idx = find(cell2mat(Water{1,3}(:,1)) == par.P);
    if any(idx)
        % Extract data for the species
        res.T_sat = Water{1,3}{idx,2};

        tab = Water{1,3}{idx,3};
        res.v = interp1(tab.T,tab.v,par.T);
        res.u = interp1(tab.T,tab.u,par.T);
        res.h = interp1(tab.T,tab.h,par.T);
        res.s = interp1(tab.T,tab.s,par.T);
    else
        error('Superheated vapor pressure not found.');
    end
    res.x=[];
end