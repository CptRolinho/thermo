classdef Thermo < handle
    %Thermo  Student-friendly thermodynamics interface (class façade)
    %
    % This class wraps the implementation in the +thermo package:
    %   - thermo.data.*   : database access
    %   - thermo.models.* : property models (ideal gas, liquid-vapor tables)
    %
    % Units (consistent with typical thermo tables used in Çengel):
    %   T     : K
    %   P     : kPa
    %   v     : m^3/kg
    %   u,h   : kJ/kg
    %   s     : kJ/(kg*K)
    %   cp,cv : kJ/(kg*K)
    %
    % Note: 1 kJ = 1 kPa*m^3, so ideal-gas relations are consistent when
    %       R is in kJ/(kg*K) and P is in kPa.

    %% Species information (read-only after construction)
    properties (SetAccess=private)
        species   (1,1) string
        model     (1,1) string
        formula   (1,1) string
        substance (1,1) string
        info      (1,1) string

        Mm (1,1) double = NaN   % kg/kmol
        R  (1,1) double = NaN   % kJ/(kg*K)
        T_c (1,1) double = NaN  % K
        P_c (1,1) double = NaN  % kPa
        V_c (1,1) double = NaN  % m^3/kmol
    end

    properties (Constant)
        % Universal gas constant
        Ru = 8.31447;         % kJ/(kmol*K)

        % Standard gravity
        g0 = 9.80665;         % m/s^2

        % Standard atmosphere
        P_std = 101.325;      % kPa
        T_std = 288.15;       % K

        % Temperature offsets
        T0_C = 273.15;        % K (0°C in Kelvin)

        % Stefan–Boltzmann constant
        sigma = 5.670374419e-8; % W/(m^2*K^4)

        % Unit conversions (convenience)
        kPa_per_bar = 100;    % kPa/bar
        kPa_per_atm = 101.325;% kPa/atm
        kPa_per_MPa = 1000;   % kPa/MPa
    end


    %% State variables and computed properties (user-readable)
    properties
        % Primary state
        T (1,1) double = NaN     % K
        P (1,1) double = NaN     % kPa
        v (1,1) double = NaN     % m^3/kg

        % Thermodynamic properties
        u (1,1) double = NaN     % kJ/kg
        h (1,1) double = NaN     % kJ/kg
        s (1,1) double = NaN     % kJ/(kg*K)

        % Ideal-gas-only
        cv (1,1) double = NaN    % kJ/(kg*K)
        cp (1,1) double = NaN    % kJ/(kg*K)
        k  (1,1) double = NaN    % -
        s_0 (1,1) double = NaN   % kJ/(kg*K)

        % Saturation / two-phase (liquid model)
        T_sat (1,1) double = NaN % K
        P_sat (1,1) double = NaN % kPa

        v_f (1,1) double = NaN
        v_g (1,1) double = NaN

        u_f (1,1) double = NaN
        u_g (1,1) double = NaN
        u_fg (1,1) double = NaN

        h_f (1,1) double = NaN
        h_g (1,1) double = NaN
        h_fg (1,1) double = NaN

        s_f (1,1) double = NaN
        s_g (1,1) double = NaN
        s_fg (1,1) double = NaN

        % Ice-related placeholders (kept for future extension)
        u_i (1,1) double = NaN
        u_ig (1,1) double = NaN
        h_i (1,1) double = NaN
        h_ig (1,1) double = NaN
        s_i (1,1) double = NaN
        s_ig (1,1) double = NaN

        x (1,1) double = NaN     % quality
        state (1,1) string = ""  % "compressed liquid" | "saturated mixture" | "superheated vapor" | etc.
    end

    properties (SetAccess=private)
        isStateSet (1,1) logical = false
    end

    %% Public API
    methods
        function th = Thermo(species, varargin)
            %Thermo Construct a Thermo object for a given species.
            %
            % th = Thermo("Water")
            % th = Thermo("Water",'T',373.15,'x',0.5)

            arguments
                species (1,1) {mustBeTextScalar}
            end
            arguments (Repeating)
                varargin
            end

            th.species = string(species);

            meta = thermo.data.speciesData(th.species);
            th.model     = string(meta.model);
            th.formula   = string(meta.formula);
            th.substance = string(meta.substance);
            th.info      = string(meta.info);

            th.Mm  = meta.Mm;
            th.R   = meta.R;
            th.T_c = meta.T_c;
            th.P_c = meta.P_c;
            th.V_c = meta.V_c;

            if ~isempty(varargin)
                th.setState(varargin{:});
            end
        end

        function setState(th, varargin)
            %setState Set/update the thermodynamic state and compute properties.
            %
            % Ideal gas: provide at least two of {T,P,v}
            % Liquid-vapor: provide at least two of {T,P,x,u,v,h,s}

            if isempty(varargin)
                error('Thermo:InvalidState','No state inputs provided.');
            end

            th.resetState();

            switch th.model
                case "ideal-gas"
                    th.applyIdealGas(varargin{:});
                case "liquid"
                    th.applyLiquid(varargin{:});
                otherwise
                    error('Thermo:UnsupportedModel','Unsupported model "%s".', th.model);
            end

            th.isStateSet = true;
        end
    end

    %% Private implementation
    methods (Access=private)
        function resetState(th)
            % Reset all state-dependent fields to defaults (NaN/empty)
            props = ["T","P","v","u","h","s","cv","cp","k","s_0", ...
                "T_sat","P_sat","v_f","v_g","u_f","u_g","u_fg", ...
                "h_f","h_g","h_fg","s_f","s_g","s_fg", ...
                "u_i","u_ig","h_i","h_ig","s_i","s_ig","x"];
            for p = props
                th.(p) = NaN;
            end
            th.state = "";
            th.isStateSet = false;
        end

        function applyIdealGas(th, varargin)
            p = inputParser;
            p.FunctionName = 'Thermo.setState (ideal-gas)';

            vPos = @(x) (isnan(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
            addParameter(p,'T',NaN,vPos);
            addParameter(p,'P',NaN,vPos);
            addParameter(p,'v',NaN,vPos);
            parse(p,varargin{:});

            T = p.Results.T; P = p.Results.P; v = p.Results.v;

            provided = ~isnan([T,P,v]);
            if sum(provided) < 2
                error('Thermo:InvalidState','Ideal gas requires at least two of {T,P,v}.');
            end

            if sum(provided) == 3
                vCalc = th.R*T/P;
                relErr = abs(v - vCalc) / max(eps, abs(vCalc));
                if relErr > 1e-6
                    error('Thermo:InconsistentState', ...
                        'Provided (T,P,v) inconsistent with ideal gas law (relErr=%.3g).', relErr);
                end
            end

            if isnan(T), T = P*v/th.R; end
            if isnan(P), P = th.R*T/v; end
            if isnan(v), v = th.R*T/P; end

            th.T = T; th.P = P; th.v = v;

            ig = thermo.models.idealGas(th.species, th.T);
            th.cv  = ig.cv;
            th.cp  = ig.cp;
            th.k   = ig.k;
            th.h   = ig.h;
            th.u   = ig.u;
            th.s_0 = ig.s_0;

            th.state = "single-phase";
        end

        function applyLiquid(th, varargin)
            % minimal validation (more detailed validation lives in thermo.models.liquid)
            p = inputParser;
            p.FunctionName = 'Thermo.setState (liquid)';

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

            vals = [p.Results.T, p.Results.P, p.Results.x, p.Results.u, p.Results.v, p.Results.h, p.Results.s];
            if sum(~isnan(vals)) < 2
                error('Thermo:InvalidState', ...
                    'Liquid model requires at least two of {T,P,x,u,v,h,s}.');
            end

            liq = thermo.models.liquid(th.species, varargin{:});

            % Core
            th.T = liq.T; th.P = liq.P; th.v = liq.v;
            th.u = liq.u; th.h = liq.h; th.s = liq.s;

            % Saturation/two-phase
            th.T_sat = emptyToNaN(liq.T_sat);
            th.P_sat = emptyToNaN(liq.P_sat);

            th.v_f  = emptyToNaN(liq.v_f);
            th.v_g  = emptyToNaN(liq.v_g);

            th.u_f  = emptyToNaN(liq.u_f);
            th.u_g  = emptyToNaN(liq.u_g);
            th.u_fg = emptyToNaN(liq.u_fg);

            th.h_f  = emptyToNaN(liq.h_f);
            th.h_g  = emptyToNaN(liq.h_g);
            th.h_fg = emptyToNaN(liq.h_fg);

            th.s_f  = emptyToNaN(liq.s_f);
            th.s_g  = emptyToNaN(liq.s_g);
            th.s_fg = emptyToNaN(liq.s_fg);

            % Ice placeholders (may be empty)
            th.u_i  = emptyToNaN(liq.u_i);
            th.u_ig = emptyToNaN(liq.u_ig);
            th.h_i  = emptyToNaN(liq.h_i);
            th.h_ig = emptyToNaN(liq.h_ig);
            th.s_i  = emptyToNaN(liq.s_i);
            th.s_ig = emptyToNaN(liq.s_ig);

            th.x     = emptyToNaN(liq.x);
            th.state = string(liq.state);
        end
    end

    %% Student-facing discovery helpers
    methods (Static)
        function constants()
            %constants Print thermodynamic/constants used across the toolkit.

            rows = {
                "Ru",         Thermo.Ru,         "kJ/(kmol*K)",  "Universal gas constant"
                "g0",         Thermo.g0,         "m/s^2",        "Standard gravity"
                "P_std",      Thermo.P_std,      "kPa",          "Standard atmospheric pressure"
                "T_std",      Thermo.T_std,      "K",            "Standard temperature"
                "T0_C",       Thermo.T0_C,       "K",            "0°C in Kelvin"
                "sigma",      Thermo.sigma,      "W/(m^2*K^4)",  "Stefan–Boltzmann constant"
                "kPa_per_bar",Thermo.kPa_per_bar,"kPa/bar",      "Pressure conversion"
                "kPa_per_atm",Thermo.kPa_per_atm,"kPa/atm",      "Pressure conversion"
                "kPa_per_MPa",Thermo.kPa_per_MPa,"kPa/MPa",      "Pressure conversion"
                };

            fprintf('\n Thermo constants:\n')
            fprintf('%-12s %-16s %-16s %s\n', 'Symbol', 'Value', 'Unit', 'Description');
            for i = 1:size(rows,1)
                fprintf('%-12s %-16.10g %-16s %s\n', rows{i,1}, rows{i,2}, rows{i,3}, rows{i,4});
            end
        end

        function substances()
            %substances Print available species from the database
            T = thermo.data.listSpecies();
            fprintf('\n Substances:\n')
            fprintf('%25s   %-10s   %-12s   %s\n','Substance','Species','Model','Information');
            for i = 1:height(T)
                fprintf('%25s | %-10s | %-12s | %s\n', ...
                    string(T.Substance{i}), string(T.Species{i}), ...
                    string(T.Model{i}), string(T.Information{i}));
            end
        end

        function methods()
            names = {'Thermo','setState','substances','properties','methods','constants'};
            text  = {'Constructor (loads species + optional state)'
                'Set/update thermodynamic state'
                'List available substances'
                'List Thermo properties (symbol, meaning, unit)'
                'List Thermo methods'
                'List Thermo constants'};
            fprintf('\n Methods:\n')
            fprintf('%12s   %-60s\n','Name','Description');
            for i = 1:numel(names)
                fprintf('%12s : %-60s\n',names{i},text{i});
            end
            fprintf('\nType: help Thermo.<method> for more help\n')
        end

        function properties()
            %properties Print curated list of properties with units/descriptions

            fprintf('%12s\n','Species Information:')
            fprintf('%10s   %-60s\n','Symbol','Description');
            specNames = {'species','model','formula','substance','info','Ru','Mm','R','T_c','P_c','V_c'};
            specText = {
                'Species identifier (key used in database)'
                'Thermodynamic model (e.g., ideal-gas, liquid)'
                'Chemical formula'
                'Common name'
                'Database info/notes'
                'Universal gas constant [kJ/(kmol*K)]'
                'Molar mass [kg/kmol]'
                'Specific gas constant [kJ/(kg*K)]'
                'Critical temperature [K]'
                'Critical pressure [kPa]'
                'Critical molar volume [m^3/kmol]'
                };
            for i = 1:numel(specNames)
                fprintf('%10s : %-60s\n', specNames{i}, specText{i});
            end

            fprintf('\n Thermodynamic state variables and properties:\n')
            fprintf('%10s   %-60s\n','Symbol','Description');
            names = {'T','P','v','u','h','s','cv','cp','k','s_0', ...
                'T_sat','P_sat','v_f','v_g','u_f','u_g','u_fg', ...
                'h_f','h_g','h_fg','s_f','s_g','s_fg','x','state'};
            text = {
                'Temperature [K]'
                'Pressure [kPa]'
                'Specific volume [m^3/kg]'
                'Specific internal energy [kJ/kg]'
                'Specific enthalpy [kJ/kg]'
                'Specific entropy [kJ/(kg*K)]'
                'Isochoric heat capacity [kJ/(kg*K)] (ideal-gas)'
                'Isobaric heat capacity [kJ/(kg*K)] (ideal-gas)'
                'Adiabatic index [-] (ideal-gas)'
                'Entropy function s0 [kJ/(kg*K)] (ideal-gas)'
                'Saturation temperature [K] (two-phase)'
                'Saturation pressure [kPa] (two-phase)'
                'Sat. liquid specific volume [m^3/kg]'
                'Sat. vapor  specific volume [m^3/kg]'
                'Sat. liquid internal energy [kJ/kg]'
                'Sat. vapor  internal energy [kJ/kg]'
                'Latent internal energy u_fg [kJ/kg]'
                'Sat. liquid enthalpy [kJ/kg]'
                'Sat. vapor  enthalpy [kJ/kg]'
                'Latent enthalpy h_fg [kJ/kg]'
                'Sat. liquid entropy [kJ/(kg*K)]'
                'Sat. vapor  entropy [kJ/(kg*K)]'
                'Latent entropy s_fg [kJ/(kg*K)]'
                'Quality x [-]'
                'State description (compressed/saturated/superheated)'
                };
            for i = 1:numel(names)
                fprintf('%10s : %-60s\n', names{i}, text{i});
            end
        end
    end
end

%% Local helper (file scope)
function y = emptyToNaN(x)
    if isempty(x)
        y = NaN;
    else
        y = double(x);
    end
end
