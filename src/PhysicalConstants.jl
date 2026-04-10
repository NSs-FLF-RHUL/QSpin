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
planck_constant = 6.62607015e-34
hbar = planck_constant / (2*pi)

# Mass of the sun; kg
mass_sun = 1.9891e30

# Speed of light in vacuum; m / s
speed_of_light_vacuum = 2.99792458e8

# Universal gravitational constant; m^3 / kg s^2
gravitational_constant = 6.67408e-11

# Electron volt; J
electron_volt = 1.60218e-19
# Electron volt in c=1 units; kg
eV_ceq1 = electron_volt / (speed_of_light_vacuum^2)

# Kiloparsec in metres; m
kiloparsec_in_m = 3.08567758e19

# Neutron mass; kg
neutron_mass = 1.674927471e-27

# Gigayear; s
# NOTE: aren't there 10^16 seconds in a giga year?
# 1 yr ~ 3e7 seconds, so 1 Gyr ~ 3e16 seconds?
gyr = 3.1556926e9

end
