from simple_tritium_transport_model import ureg, Model
import numpy as np
from helpers import (
    substract_background_from_measurements,
    cumulative_activity,
    background_sub,
)


background_1 = 0.282 * ureg.Bq
background_2 = 0.272 * ureg.Bq
raw_measurements = {
    1: {
        1: 1.102 * ureg.Bq,
        2: 0.300 * ureg.Bq,
        3: 0.347 * ureg.Bq,
        4: 0.308 * ureg.Bq,
        "background": background_1,
    },
    2: {
        1: 1.412 * ureg.Bq,
        2: 0.352 * ureg.Bq,
        3: 0.354 * ureg.Bq,
        4: 0.329 * ureg.Bq,
        "background": background_1,
    },
    3: {
        1: 1.069 * ureg.Bq,
        2: 0.447 * ureg.Bq,
        3: 0.662 * ureg.Bq,
        4: 0.347 * ureg.Bq,
        "background": background_1,
    },
    4: {
        1: 2.415 * ureg.Bq,
        2: 0.569 * ureg.Bq,
        3: 0.570 * ureg.Bq,
        4: 0.333 * ureg.Bq,
        "background": background_1,
    },
    5: {
        1: 2.538 * ureg.Bq,
        2: 0.510 * ureg.Bq,
        3: 0.348 * ureg.Bq,
        4: 0.292 * ureg.Bq,
        "background": background_2,
    },
    6: {
        1: 1.188 * ureg.Bq,
        2: 0.351 * ureg.Bq,
        3: 0.343 * ureg.Bq,
        4: 0.313 * ureg.Bq,
        "background": background_2,
    },
    7: {
        1: 1.163 * ureg.Bq,
        2: 0.541 * ureg.Bq,
        3: 0.349 * ureg.Bq,
        4: 0.316 * ureg.Bq,
        "background": background_2,
    },
}

measurements_after_background_sub = substract_background_from_measurements(
    raw_measurements
)


# time starts at 04/03 10:20 AM
# 04/04 10:20 AM = 24 hours = 1 * ureg.day + 0 * ureg.hour + 0 * ureg.minute
# 04/05 10:20 AM = 48 hours = 2 * ureg.day + 0 * ureg.hour + 0 * ureg.minute
replacement_times = [
    # 04/03 22:27
    0 * ureg.day + 12 * ureg.hour + 7 * ureg.minute,
    # 04/04 10:06
    1 * ureg.day + 0 * ureg.hour - 14 * ureg.minute,
    # 04/04 23:39
    1 * ureg.day + 13 * ureg.hour + 19 * ureg.minute,
    # 04/05 15:15
    2 * ureg.day + 4 * ureg.hour + 55 * ureg.minute,
    # 04/06 17:46
    3 * ureg.day + 7 * ureg.hour + 26 * ureg.minute,
    # 04/07 15:33
    4 * ureg.day + 5 * ureg.hour + 13 * ureg.minute,
    # 04/09 11:34
    6 * ureg.day + 1 * ureg.hour + 14 * ureg.minute,
]

replacement_times = sorted(replacement_times)

# # Cumulative values

cumulative_release = cumulative_activity(measurements_after_background_sub)

# Model

baby_diameter = 1.77 * ureg.inches - 2 * 0.06 * ureg.inches  # from CAD drawings
baby_radius = 0.5 * baby_diameter
baby_volume = 0.125 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section
calculated_TBR = 4.57e-4 * ureg.particle * ureg.neutron**-1  # stefano 1/22/2024
baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=calculated_TBR,
)


mass_transport_coeff_factor = 3

baby_model.k_top *= mass_transport_coeff_factor * 0.09
optimised_ratio = 0.29
baby_model.k_wall = baby_model.k_top * optimised_ratio

exposure_time = 12 * ureg.hour

baby_model.irradiations = [
    [0 * ureg.hour, 0 + exposure_time],
    [24 * ureg.hour, 24 * ureg.hour + exposure_time],
]

# calculated from Kevin's activation foil analysis
P383_neutron_rate = 4.95e8 * ureg.neutron * ureg.s**-1
A325_neutron_rate = 2.13e8 * ureg.neutron * ureg.s**-1

neutron_rate_relative_uncertainty = 0.089
baby_model.neutron_rate = P383_neutron_rate + A325_neutron_rate
