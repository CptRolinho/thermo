function par = idealGas(species, T)
%thermo.models.idealGas  Ideal gas property evaluation using thermo.mat payload
% Parameters for thermodynamic model for a pseudo species as an ideal gas
% Reference: Çengel
% v.0.0.1
%
% par = thermo.models.idealGas("Air", 300)

    arguments
        species (1,1) {mustBeTextScalar}
        T (1,1) double {mustBePositive, mustBeFinite}
    end

    payload = thermo.data.getPayload(species); % e.g., Air cell
    if ~iscell(payload) || numel(payload) < 2
        error('thermo:models:InvalidPayload', ...
            'Ideal gas payload for "%s" must be a cell with at least 2 elements.', species);
    end

    tab1 = payload{1};
    tab2 = payload{2};

    % Use safe interpolation with bounds error (recommended)
    par.cv  = interp1safe(tab1.T, tab1.cv,  T);
    par.cp  = interp1safe(tab1.T, tab1.cp,  T);
    par.k   = interp1safe(tab1.T, tab1.k,   T);
    par.h   = interp1safe(tab2.T, tab2.h,   T);
    par.u   = interp1safe(tab2.T, tab2.u,   T);
    par.s_0 = interp1safe(tab2.T, tab2.s_0, T);
end
