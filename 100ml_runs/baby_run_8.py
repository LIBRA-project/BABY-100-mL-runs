from libra_toolbox.tritium.model import ureg, Model
import numpy as np
from libra_toolbox.tritium.helpers import (
    substract_background_from_measurements,
    cumulative_activity,
    background_sub,
)


BL_1 = 0.334 * ureg.Bq
BL_2 = 0.322 * ureg.Bq
background = (BL_1 + BL_2) / 2

raw_measurements = {
    1: {
        1: 0.396 * ureg.Bq,
        2: 0.375 * ureg.Bq,
        3: 4.575 * ureg.Bq,
        4: 0.448 * ureg.Bq,
        "background": background,
    },
    2: {
        1: 0.386 * ureg.Bq,
        2: 0.417 * ureg.Bq,
        3: 5.659 * ureg.Bq,
        4: 0.509 * ureg.Bq,
        "background": background,
    },
    3: {
        1: 0.393 * ureg.Bq,
        2: 0.410 * ureg.Bq,
        3: 6.811 * ureg.Bq,
        4: 0.492 * ureg.Bq,
        "background": background,
    },
    4: {
        1: 0.406 * ureg.Bq,
        2: 0.403 * ureg.Bq,
        3: 4.864 * ureg.Bq,
        4: 0.467 * ureg.Bq,
        "background": background,
    },
    5: {
        1: 0.322 * ureg.Bq,
        2: 0.369 * ureg.Bq,
        3: 1.900 * ureg.Bq,
        4: 0.470 * ureg.Bq,
        "background": background,
    },
    6: {
        1: 0.343 * ureg.Bq,
        2: 0.363 * ureg.Bq,
        3: 0.492 * ureg.Bq,
        4: 0.361 * ureg.Bq,
        "background": background,
    },
    7: {
        1: 0.287 * ureg.Bq,
        2: 0.367 * ureg.Bq,
        3: 0.353 * ureg.Bq,
        4: 0.328 * ureg.Bq,
        "background": background,
    },
}

measurements_after_background_sub = substract_background_from_measurements(
    raw_measurements
)


# time starts at 7/29/2024 9:28:00 AM
replacement_times = [
    0.51 * ureg.day,
    1.00 * ureg.day,
    1.52 * ureg.day,
    2.02 * ureg.day,
    3.10 * ureg.day,
    4.08 * ureg.day,
    6.24 * ureg.day,
    # 9.15 * ureg.day,
]

replacement_times = sorted(replacement_times)

# # Cumulative values

cumulative_release = cumulative_activity(measurements_after_background_sub)


# Model

baby_diameter = 1.77 * ureg.inches - 2 * 0.06 * ureg.inches  # from CAD drawings
baby_radius = 0.5 * baby_diameter
baby_volume = 0.100 * ureg.L  # TODO double check this value
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section

# This is just adapted for fitting, not an actual OpenMC calculation yet
calculated_TBR = (
    5.4e-4 * ureg.particle * ureg.neutron**-1
)  # TODO double check this value with Stefano


mass_transport_coeff_factor = 3

k_top = 4.9e-7 * ureg.m * ureg.s**-1 * mass_transport_coeff_factor * 1.15
optimised_ratio = 3e-2
k_wall = k_top * optimised_ratio

exposure_time = 12 * ureg.hour

irradiations = [
    [0 * ureg.hour, 0 + exposure_time],
    [24 * ureg.hour, 24 * ureg.hour + exposure_time],
]

# calculated from Kevin's activation foil analysis
P383_neutron_rate = 4.95e8 * ureg.neutron * ureg.s**-1
A325_neutron_rate = 2.13e8 * ureg.neutron * ureg.s**-1

neutron_rate_relative_uncertainty = 0.089
neutron_rate = (
    P383_neutron_rate + A325_neutron_rate
) / 2  # the neutron rate is divided by two to acount for the double counting (two detectors)

baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=calculated_TBR,
    k_top=k_top,
    k_wall=k_wall,
    neutron_rate=neutron_rate,
    irradiations=irradiations,
)
