"""
Defines a number of fundamental physical constants for use in simulations.
"""
module PhysicalConstants
using PhysicalConstants.CODATA2022:
    NewtonianConstantOfGravitation,
    SpeedOfLightInVacuum,
    PlanckConstant,
    ReducedPlanckConstant,
    NeutronMass,
    ElementaryCharge
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
const planck_constant = PlanckConstant.val
const hbar = ReducedPlanckConstant.val

# Mass of the sun; kg
const mass_sun = 1.98855e30

# Speed of light in vacuum; m / s
const speed_of_light_vacuum = SpeedOfLightInVacuum.val

# Universal gravitational constant; m^3 / kg s^2
const gravitational_constant = NewtonianConstantOfGravitation.val

# Electron volt; J (numerically same as elementary charge in Coulomb)
const electron_volt = ElementaryCharge.val
# Electron volt in c=1 units; kg
const eV_ceq1 = electron_volt / (speed_of_light_vacuum^2)

# Kiloparsec in metres; m
const kiloparsec_in_m = 3.08567758e19

# Neutron mass; kg
const neutron_mass = NeutronMass.val
const atomic_mass = AtomicMassConstant.val 

# Gigayear; s
const gyr = 86400 * 365.24219 * 1e9

end
