# thermo
MATLAB class for managing and calculating thermodynamic properties of substances

## Overview
The `thermo` class is designed for managing and calculating thermodynamic properties of substances. It provides simple and robust tools for calculating and retrieving various thermodynamic properties based on the thermodynamics property tables from Çengel, Y. A., Boles, M. A., & Kanoğlu, M. (2024). Thermodynamics: An engineering approach, using specified thermodynamic models and state variables such as temperature (T), pressure (P), and specific volume (v).

This class supports multiple species, each with data that can be fetched and interpolated based on input conditions. It is capable of handling ideal gases and liquids and can be extended to other phases or models as needed.

## Features
- **Calculation of Thermodynamic Properties**: Dynamically calculate properties based on temperature, pressure, and volume.
- **Support for Multiple Species**: Includes predefined data for common substances like Air and Water, and can be extended to include more.
- **Extensible**: Designed to be easily extended to handle other phases or thermodynamic models.
- **Direct Access**: Properties are accessible directly through object properties after initialization.

## Usage

To use the `thermo` class, you first need to initialize an object for a specific species. Here are some examples:

```matlab
% Initialize the thermo object for Air
th = thermo('Air');

% Initialize with specific state variables
th = thermo('Air', 'T', 300, 'P', 101.325);
```
## Getting Started
To get a complete list of available substances, properties, and methods, you can use the following commands within MATLAB:
```matlab
thermo.substances
thermo.properties
thermo.methods
```

## Aim
The class aims to provide accurate results consistent with known thermodynamic models and literature values, making it suitable for both engineering and scientific applications.

## Contributions
Contributions to extend the functionality, add new species, or improve existing models are welcome. Please feel free to fork the repository, make your changes, and submit a pull request.
