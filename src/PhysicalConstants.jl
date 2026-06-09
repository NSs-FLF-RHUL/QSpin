"""
Defines a number of fundamental physical constants for use in simulations.
"""
module PhysicalConstants

export planck_constant,
    hbar,
    mass_sun,
    speed_of_light_vacuum,
    gravitational_constant,
    electron_volt,
    eV_ceq1,
    kiloparsec_in_m,
    neutron_mass,
    gyr

# planck constant; m^2 kg / s
const planck_constant = 6.62607015e-34
const hbar = planck_constant / (2*pi)

# Mass of the sun; kg
const mass_sun = 1.98855e30

# Speed of light in vacuum; m / s
const speed_of_light_vacuum = 299792458.0

# Universal gravitational constant; m^3 / kg s^2
const gravitational_constant = 6.67408e-11

# Electron volt; J
const electron_volt = 1.66053906892e-19
# Electron volt in c=1 units; kg
const eV_ceq1 = electron_volt / (speed_of_light_vacuum^2)

# Kiloparsec in metres; m
const kiloparsec_in_m = 3.08567758e19

# Neutron mass; kg
const neutron_mass = 1.66053906892e-27

# Gigayear; s
# NOTE: aren't there 10^16 seconds in a giga year?
# 1 yr ~ 3e7 seconds, so 1 Gyr ~ 3e16 seconds?
const gyr = 3.1556926e9


# Unit Converters for SI and CGS units for commonly used NS parameters in the literature

const kg2g = 1e3
const g2kg = 1e-3
const cm2m = 1e-2
const m2cm = 1e2
const fm2m = 1e-15
const m2fm = 1e15
const dyne2N = 1e-5
const N2dyne = 1e5

end
