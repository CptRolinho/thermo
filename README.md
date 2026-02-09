# Thermo — MATLAB Thermodynamics Toolkit

A lightweight, research-grade **thermodynamics toolkit for MATLAB**, designed to be **easy for students** and **powerful for researchers**.

Thermo is based on thermodynamic property tables from:
> Çengel, Y. A., Boles, M. A., & Kanoğlu, M. (2024). *Thermodynamics: An Engineering Approach*.

It combines:
- a **user-friendly class API** (`Thermo`) for interactive use, and
- a **function-based package API** (`+thermo`) for fast loops, sweeps, and research workflows.

The physics and numerical logic live in the package; the class is a clean façade on top.

---

## Key Features

- Ideal-gas and liquid–vapor property evaluation (Air, Water, R134a)
- Automatic state resolution (compressed / saturated / superheated)
- Robust interpolation with protection against incomplete tables
- Consistent SI engineering units (kPa, K, kJ/kg, m³/kg)
- Works seamlessly in scripts, live scripts, and large parametric loops
- Extensible architecture for new fluids and models

---

## Folder Structure

```
project_root/
│
├─ Thermo.m                % Main class (user-facing API)
├─ setupThermo.m           % One-command setup & verification
│
├─ +thermo/                % Package (implementation layer)
│   ├─ +data/              % Databases and metadata
│   │   ├─ thermo.mat
│   │   ├─ speciesData.m
│   │   ├─ getPayload.m
│   │   └─ listSpecies.m
│   │
│   └─ +models/            % Thermodynamic models
│       ├─ idealGas.m
│       ├─ liquid.m
│       └─ private/
│           └─ interp1safe.m
│
├─ tests/
│   └─ smokeTest_Thermo.m
│
└─ README.md
```

---

## Installation & Setup (Recommended)

From the project root, run once per MATLAB session:

```matlab
setupThermo
```

This will:
- add the correct folders to the MATLAB path,
- verify that the class and package are available,
- run a smoke test to confirm everything works.

---

## Quick Start

### Create a thermodynamic state

```matlab
th = Thermo("Water", 'T', 373.15, 'x', 0.5);
```

### Inspect properties

```matlab
th.state    % "saturated mixture"
th.P        % saturation pressure [kPa]
th.h        % enthalpy [kJ/kg]
th.s        % entropy [kJ/(kg·K)]
```

### Update the state

```matlab
th.setState('P', 200, 'T', 450);
```

All properties are updated automatically.

---

## Supported Input Combinations

### Ideal Gas (e.g. Air)
Provide **any two** of:
- `T` (K)
- `P` (kPa)
- `v` (m³/kg)

Example:
```matlab
air = Thermo("Air", 'T', 300, 'P', 101.325);
```

### Liquid–Vapor Fluids (Water, R134a)
Provide **any two** of:
- `T` (K)
- `P` (kPa)
- `x` (quality)
- `u` (kJ/kg)
- `v` (m³/kg)
- `h` (kJ/kg)
- `s` (kJ/(kg·K))

The thermodynamic state (compressed / saturated / superheated) is inferred automatically.

---

## Research / Loop-Based Usage

For parameter sweeps, optimization, or post-processing, call the **package models directly**:

```matlab
par = thermo.models.liquid("Water", 'T', 373.15, 'x', 0.5);
```

Example loop:
```matlab
Ts = linspace(300, 500, 50);
hs = zeros(size(Ts));

for i = 1:numel(Ts)
    p = thermo.models.liquid("Water", 'T', Ts(i), 'P', 200);
    hs(i) = p.h;
end
```

Returned values are plain MATLAB structs—easy to store, vectorize, and analyze.

---

## Discoverability & Help

```matlab
Thermo.substances()   % list available fluids
Thermo.methods()      % list class methods
Thermo.properties()   % list properties with units
Thermo.constants()    % list thermodynamic constants
```

---

## Units

All calculations use a consistent SI engineering unit system:

| Quantity | Unit |
|--------|------|
| Temperature | K |
| Pressure | kPa |
| Specific volume | m³/kg |
| Enthalpy | kJ/kg |
| Internal energy | kJ/kg |
| Entropy | kJ/(kg·K) |
| Gas constant | kJ/(kg·K) |

Note: `1 kJ = 1 kPa·m³`, so ideal-gas relations are internally consistent.

---

## Design Philosophy

- **Single source of truth**: all physics lives in `+thermo/models`
- **Thin class façade**: `Thermo` only manages state and validation
- **Safe numerics**: no silent extrapolation; incomplete tables are handled robustly
- **Pedagogical first**: clear errors, discoverable API, textbook consistency

---

## Extending the Toolkit

- Add a new fluid → update `thermo.mat` + `speciesData.m`
- Add a new model → create `thermo.models.<model>.m`
- Add plots → new methods in `Thermo` (e.g. `plotTS`, `plotPH`)

---

## Disclaimer

This code is intended for **educational and research use**.
Thermodynamic data are sourced from standard textbook tables and should be validated against authoritative references for critical applications.

---

Happy thermodynamics 🚀

