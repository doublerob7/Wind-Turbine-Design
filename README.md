#Wind Turbine Design

## Description
This repo contains the Python 2.7 code developed for ME 4470 Wind and Ocean Energy, taught Fall 2015. The code handles the design and analysis of a wind turbine plant through several intermediary steps starting with wind and airfoil data.

###Step 1 - Reading and formatting raw wind data
Real, logged wind data recorded at Laramie Regional Airport is read and formatted for use.
###Step 2 - Statistical analysis of wind data
Wind data is analyzed over specific intervals to characterize wind conditions at site. Probability Density Functions are calculated to determine mean wind speed and variance. Weibull and Rayleigh distributions are fit to the data.
###Step 3 - Blade Element Momentum Analysis
A wind turbine blade is simulated as a series of blade elements. Each element undergoes an aerodynamic load analysis based on incoming wind velocity and aerodynamic data for the section profile to determine the forces acting on the section. Turbine Thrust, Torque and Power are then calculated based on the blade forces from 3 blades.
###Step 4 - Blade Shear/Moment and Deformation
The blades are treated as I-beams to utilize beam theory to determine shear forces, moments and deflections of the blades.
###Step 5 - Turbine/Generator system Power Curve
The turbine hub (3 blades) is coupled to a generator, with given properties, and this system is analyzed to produce a Power Curve. The Power Curve describes the electrical power produced as a function of incoming wind velocity.
###Step 6 - Blade Pitch Control Schedule
The pitch of the turbine blades is considered in order to broaden the range of power-producing wind speeds. Blade pitch is incremented and a Power Curve is produced for each pitch value; each pitch corresponding to maximum generator output power is calculated. A Pitch Schedule function is then fit to these values and describes the pitch to be set by a controller for maximum power output over the entire range of operatational windspeeds.
###Step 7 - Turbine Plant and Storage based on simulated load
The turbine model is now considered in the context of a wind plant, consisting of several turbines to meet the power demands of a simulated city. The city's power consumption is modeled as a full-high full-low cycle over the course of the day; full-high for daylight hours, full-low for night hours. A storage capacity is considered to smooth out the highly variable wind power, and a minimum storage capacity is calculated to provide continuous power to the city.

## Roadmap - TODO
- [ ] Complete step 2, statistical analysis
- [ ] Complete step 4, blade loads
- [ ] Model wind turbine in CAD