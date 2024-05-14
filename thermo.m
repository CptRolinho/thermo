classdef    thermo < handle
    %THERMO: Class for managing and calculating thermodynamic properties of substances
    % This class provides tools for calculating and retrieving various
    % thermodynamic properties based on the specified thermodynamic model
    % and state variables such as temperature (T), pressure (P), and specific volume (v).
    % It supports multiple species, each with data that can be fetched and
    % interpolated based on the input conditions.
    %
    % The class handles ideal gases and liquids and can be extended to other
    % phases or models as needed. The properties are accessible directly
    % through object properties after initialization and setting state conditions.
    %
    % Usage:
    %   th = thermo(species); % Initializes the object for a given species.
    %   th = thermo(species, 'T', 300, 'P', 101.325); % Example of specifying state.
    %
    % Supported species include but are not limited to Air and Water. Each species
    % must have predefined data that includes temperature-dependent properties
    % such as cp (specific heat at constant pressure), cv (specific heat at constant volume),
    % and other thermodynamic properties.
    %
    % To see a complete list of substances, properties and methods:
    %   thermo.substances
    %   thermo.properties
    %   thermo.methods
    %
    % The class aims to provide accurate results consistent with known thermodynamic
    % models and literature values, suitable for engineering and scientific applications.
    %
    % Created by: Miguel Rolinho Clemente, May 2024

    properties      % Species information
        species;    % Identifier of current species (mix of name, formula, or more common identifier)
        model;      % Thermodynamic model (e.g., 'ideal-gas', 'liquid')
        formula;    % Chemical formula of the species (e.g., 'H2O', 'N2')
        substance;  % Common name of the species (e.g., 'Air')
        Ru = 8.31447; % Universal gas constant (kJ/kmol*K)
        Mm;         % Molar mass (kg/kmol)
        R;          % Individual gas constant (kJ/kg*K)
        T_c;        % Critical-point temperature (K)
        P_c;        % Critical-point pressure
        V_c;        % Critical-point molar volume (m3/kmol)
    end

    properties      % Thermodynamic state variables and properties
        T;          % Temperature (K)
        P;          % Pressure (kPa)
        v;          % Specific volume (m3/kg)
        u;          % Internal energy (kJ/kg)
        h;          % Enthalpy (kJ/kg)
        s;          % Entropy (kJ/kg*K)
        cv;         % Specific heat capacity at constant volume (kJ/kg*K)
        cp;         % Specific heat capacity at constant pressure (kJ/kg*K)
        k;          % Adiabatic index
        s_0;        % Entropy function zero (kJ/kg*K)
        T_sat;      % Saturated temperature (K)
        P_sat;      % Saturated pressure (kPa)
        v_f;        % Saturated liquid specific volume (m3/kg)
        v_g;        % Saturated vapor specific volume (m3/kg)
        u_f;        % Saturated liquid internal energy (kJ/kg)
        u_g;        % Saturated vapor internal energy (kJ/kg)
        u_fg;       % Evaporation internal energy (kJ/kg)
        u_i;        % Saturated ice internal energy (kJ/kg)
        u_ig;       % Sublimation internal energy (kJ/kg)
        h_f;        % Saturated liquid enthalpy (kJ/kg)
        h_g;        % Saturated vapor enthalpy (kJ/kg)
        h_fg;       % Evaporation enthalpy (kJ/kg)
        h_i;        % Saturated ice enthalpy (kJ/kg)
        h_ig;       % Sublimation enthalpy (kJ/kg)
        s_f;        % Saturated liquid entropy (kJ/kg*K)
        s_g;        % Saturated vapor entropy (kJ/kg*K)
        s_fg;       % Evaporation entropy (kJ/kg*K)
        s_i;        % Saturated ice entropy (kJ/kg*K)
        s_ig;       % Sublimation entropy (kJ/kg*K)
        x;          % Quality, fraction of vapor in a saturated mixture
        state;      % State of liquid-vapor mixture

    end
    methods
        function th = thermo(species, varargin)
            % THERMO: Thermodynamic object constructor for the thermo class
            % Initializes the thermodynamic properties for a given species.
            % Usage: th = thermo(species, varargin);
            % Inputs:
            %   species: String, specifies the chemical species
            %   varargin: Name-value pairs for setting initial state conditions
            %       T, P, v, u, x

            % Load species specific data from an external function
            th.species = species;
            par = speciesData(th.species);

            % Assign general species properties
            th.model = par.model;
            th.formula = par.formula;
            th.substance = par.substance;
            th.Mm = par.Mm;     % Molar mass (kg/kmol)
            th.R = par.R;       % Gas constant (kJ/kg*K)
            th.T_c = par.T_c;   % Critical-point temperature (K)
            th.P_c = par.P_c;   % Critical-point pressure
            th.V_c = par.V_c;   % Critical-point molar volume (m3/kmol)

            % Initialize and parse additional input parameters
            switch th.model
                case 'ideal-gas'
                    p = inputParser;
                    addParameter(p, 'T', NaN);
                    addParameter(p, 'P', NaN);
                    addParameter(p, 'v', NaN);
                    parse(p, varargin{:});

                    th.T = p.Results.T;
                    th.P = p.Results.P;
                    th.v = p.Results.v;

                    % Ensure sufficient parameters are provided to define state
                    provided = ~isnan([th.T, th.P, th.v]);
                    paramsProvided = sum(provided);
                    if paramsProvided < 2
                        error('Insufficient inputs: Please provide at least two variables among T, P, v.');
                    end

                    % Calculate properties based on ideal gas law
                    if provided(1) && provided(2) % T and P provided
                        th.v = th.R*th.T/th.P;
                    elseif provided(1) && provided(3) % T and v provided
                        th.P = th.R*th.T/th.v;
                    elseif provided(2) && provided(3) % P and v provided
                        th.T = th.P*th.v/th.R;
                    end

                    % Load temperature-dependent properties
                    par = ideal_gas(species,th.T);
                    th.cv = par.cv;     % Specific heat capacity at constant volume (kJ/kg*K)
                    th.cp = par.cp;     % Specific heat capacity at constant pressure (kJ/kg*K)
                    th.k = par.k;       % Adiabatic index
                    th.h = par.h;       % Enthalpy (kJ/kg)
                    th.u = par.u;       % Internal energy (kJ/kg)
                    th.s_0 = par.s_0;   % Entropy function zero (kJ/kg*K)

                case 'liquid'
                    % Initialize and parse additional input parameters
                    p = inputParser;
                    addParameter(p, 'T', NaN);
                    addParameter(p, 'P', NaN);
                    addParameter(p, 'x', NaN);
                    addParameter(p, 'u', NaN);
                    parse(p, varargin{:});

                    th.T = p.Results.T;
                    th.P = p.Results.P;
                    th.x = p.Results.x;
                    th.u = p.Results.u;

                    % Ensure sufficient parameters are provided to define state
                    provided = ~isnan([th.T, th.P, th.x, th.u]);
                    paramsProvided = sum(provided);
                    if paramsProvided < 2
                        error('Insufficient inputs: Please provide at least two variables among T, P, x, u.');
                    end
                    par = feval(['liquid_',species],'T',th.T,'P',th.P,'x',th.x,'u',th.u);

                    % % Load temperature, or pressure, dependent properties
                    % if provided(1) && provided(2) % T and P provided
                    %     par = feval(['liquid_',species],'T',th.T,'P',th.P);
                    % elseif provided(1) && provided(3) % T and x provided
                    %     par = feval(['liquid_',species],'T',th.T,'x',th.x);
                    % elseif provided(1) && provided(4) % T and u provided
                    %     par = feval(['liquid_',species],'T',th.T,'u',th.u);
                    % elseif provided(2) && provided(3) % P and x provided
                    %     par = feval(['liquid_',species],'P',th.P,'x',th.x);
                    % elseif provided(2) && provided(4) % P and u provided
                    %     par = feval(['liquid_',species],'P',th.P,'P',th.u);
                    % end

                    th.T     = par.T;        % Temperature (K)
                    th.P     = par.P;        % Pressure (kPa)
                    th.v     = par.v;        % Specific volume (m3/kg)
                    th.u     = par.u;        % Internal energy (kJ/kg)
                    th.h     = par.h;        % Enthalpy (kJ/kg)
                    th.s     = par.s;        % Entropy (kJ/kg*K)
                    th.T_sat = par.T_sat;    % Saturated temperature (K)
                    th.P_sat = par.P_sat;    % Saturated pressure (kPa)
                    th.v_f   = par.v_f;      % Saturated liquid specific volume (m3/kg)
                    th.v_g   = par.v_g;      % Saturated vapor specific volume (m3/kg)
                    th.u_f   = par.u_f;      % Saturated liquid internal energy (kJ/kg)
                    th.u_g   = par.u_g;      % Saturated vapor internal energy (kJ/kg)
                    th.u_fg  = par.u_fg;     % Evaporation internal energy (kJ/kg)
                    th.h_f   = par.h_f;      % Saturated liquid enthalpy (kJ/kg)
                    th.h_g   = par.h_g;      % Saturated vapor enthalpy (kJ/kg)
                    th.h_fg  = par.h_fg;     % Evaporation enthalpy (kJ/kg)
                    th.s_f   = par.s_f;      % Saturated liquid entropy (kJ/kg*K)
                    th.s_g   = par.s_g;      % Saturated vapor entropy (kJ/kg*K)
                    th.s_fg  = par.s_fg;     % Evaporation entropy (kJ/kg*K)
                    th.u_i   = par.u_i;      % Saturated ice internal energy (kJ/kg)
                    th.u_ig  = par.u_ig;     % Sublimation internal energy (kJ/kg)
                    th.h_i   = par.h_i;      % Saturated ice enthalpy (kJ/kg)
                    th.h_ig  = par.h_ig;     % Sublimation enthalpy (kJ/kg)
                    th.s_i   = par.s_i;      % Saturated ice entropy (kJ/kg*K)
                    th.s_ig  = par.s_ig;     % Sublimation entropy (kJ/kg*K)
                    th.x     = par.x;        % Quality
                    th.state = par.state;

                case 'none'
                otherwise
                    error('Model not found or not supported.')
            end
        end
    end
    methods(Static)
        function properties
            % PROPERTIES: Displays a list of all propetries of the thermo class
            names = {
                'species'
                'model'
                'formula'
                'substance'
                'Ru'
                'Mm'
                'R'
                'T_c'
                'P_c'
                'V_c'
                };
            text = {
                'Identifier of current species (mix of name, formula, or more common identifier)'
                'Thermodynamic model (e.g., ''ideal-gas'', ''liquid'')'
                'Chemical formula of the species (e.g., ''H2O'', ''N2'')'
                'Common name of the species (e.g., ''Air'')'
                'Universal gas constant (kJ/kmol*K)'
                'Molar mass (kg/kmol)'
                'Individual gas constant (kJ/kg*K)'
                'Critical-point temperature (K)'
                'Critical-point pressure'
                'Critical-point molar volume (m3/kmol)'
                };

            fprintf('%12s\n','Species Information:')
            fprintf('%10s   %-50s\n','Symbol','Description');
            for i = 1:length(names)
                fprintf('%10s : %-50s\n',names{i},text{i});
            end

            names = {
                'T'
                'P'
                'v'
                'u'
                'h'
                's'
                'cv'
                'cp'
                'k'
                's_0'
                'T_sat'
                'P_sat'
                'v_f'
                'v_g'
                'u_f'
                'u_g'
                'u_fg'
                'u_i'
                'u_ig'
                'h_f'
                'h_g'
                'h_fg'
                'h_i'
                'h_ig'
                's_f'
                's_g'
                's_fg'
                's_i'
                's_ig'
                'x'
                'state'
                };
            text = {
                'Temperature (K)'
                'Pressure (kPa)'
                'Specific volume (m3/kg)'
                'Internal energy (kJ/kg)'
                'Enthalpy (kJ/kg)'
                'Entropy (kJ/kg*K)'
                'Specific heat capacity at constant volume (kJ/kg*K)'
                'Specific heat capacity at constant pressure (kJ/kg*K)'
                'Adiabatic index'
                'Entropy function zero (kJ/kg*K)'
                'Saturated temperature (K)'
                'Saturated pressure (kPa)'
                'Saturated liquid specific volume (m3/kg)'
                'Saturated vapor specific volume (m3/kg)'
                'Saturated liquid internal energy (kJ/kg)'
                'Saturated vapor internal energy (kJ/kg)'
                'Evaporation internal energy (kJ/kg)'
                'Saturated ice internal energy (kJ/kg)'
                'Sublimation internal energy (kJ/kg)'
                'Saturated liquid enthalpy (kJ/kg)'
                'Saturated vapor enthalpy (kJ/kg)'
                'Evaporation enthalpy (kJ/kg)'
                'Saturated ice enthalpy (kJ/kg)'
                'Sublimation enthalpy (kJ/kg)'
                'Saturated liquid entropy (kJ/kg*K)'
                'Saturated vapor entropy (kJ/kg*K)'
                'Evaporation entropy (kJ/kg*K)'
                'Saturated ice entropy (kJ/kg*K)'
                'Sublimation entropy (kJ/kg*K)'
                'Quality, fraction of vapor in a saturated mixture'
                'State of liquid-vapor mixture'
                };
            fprintf('\n Thermodynamic state variables and properties:\n')
            fprintf('%10s   %-50s\n','Symbol','Description');
            for i = 1:length(names)
                fprintf('%10s : %-50s\n',names{i},text{i});
            end
        end

        function methods
            % METHODS: Displays a list of methods of the thermo class
            names = {'thermo','substances','properties','methods'};
            text = {'Creator function for thermodynamic object'
                'List available substances'
                'List thermo properties'
                'List thermo methods'
                };
            fprintf('\n Methods:\n')
            fprintf('%10s   %-50s\n','Name','Description');
            for i = 1:length(names)
                fprintf('%10s : %-50s\n',names{i},text{i});
            end
            fprintf('\nType: help thermo.<method> for more help\n')
        end

        function substances
            % SUBSTANCES: Displays the list of available species of the thermo class
            load('thermo.mat', 'SpeciesData');
            fprintf('\n Substances:\n')
            fprintf('%25s   %-7s   %-9s   %s\n','Substance','Species','Model','Information');
            for i = 1:height(SpeciesData)
                fprintf('%25s | %-7s | %-9s | %s\n',SpeciesData.Substance{i},SpeciesData.Species{i},SpeciesData.Model{i},SpeciesData.Information{i});
            end
        end
    end
end