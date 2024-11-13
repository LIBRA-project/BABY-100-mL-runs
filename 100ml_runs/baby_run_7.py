from libra_toolbox.tritium.model import ureg, Model
import numpy as np
from libra_toolbox.tritium.helpers import (
    substract_background_from_measurements,
    cumulative_activity,
    background_sub,
)


background_1 = 0.278 * ureg.Bq  # unclear what the background is on COUNT 2
background_2 = 0.278 * ureg.Bq
raw_measurements = {
    1: {
        1: 0.990 * ureg.Bq,
        2: 0.302 * ureg.Bq,
        3: 0.625 * ureg.Bq,
        4: 0.265 * ureg.Bq,
        "background": background_1,
    },
    2: {
        1: 1.410 * ureg.Bq,
        2: 0.337 * ureg.Bq,
        3: 0.326 * ureg.Bq,
        4: 0.303 * ureg.Bq,
        "background": background_1,
    },
    3: {
        1: 1.037 * ureg.Bq,
        2: 0.310 * ureg.Bq,
        3: 0.617 * ureg.Bq,
        4: 0.287 * ureg.Bq,
        "background": background_1,
    },
    4: {
        1: 2.457 * ureg.Bq,
        2: 0.365 * ureg.Bq,
        3: 0.429 * ureg.Bq,
        4: 0.285 * ureg.Bq,
        "background": background_1,
    },
    5: {
        1: 2.576 * ureg.Bq,
        2: 0.496 * ureg.Bq,
        3: 0.370 * ureg.Bq,
        4: 0.287 * ureg.Bq,
        "background": background_2,
    },
    6: {
        1: 1.217 * ureg.Bq,
        2: 0.365 * ureg.Bq,
        3: 0.323 * ureg.Bq,
        4: 0.282 * ureg.Bq,
        "background": background_2,
    },
    7: {
        1: 1.147 * ureg.Bq,
        2: 0.428 * ureg.Bq,
        3: 0.370 * ureg.Bq,
        4: 0.290 * ureg.Bq,
        "background": background_2,
    },
    8: {
        1: 0.940 * ureg.Bq,
        2: 0.526 * ureg.Bq,
        3: 0.422 * ureg.Bq,
        4: 0.300 * ureg.Bq,
        "background": background_2,
    },
    9: {
        1: 0.783 * ureg.Bq,
        2: 0.462 * ureg.Bq,
        3: 0.439 * ureg.Bq,
        4: 0.340 * ureg.Bq,
        "background": 0.249 * ureg.Bq,
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
    # 04/12 15:03
    9 * ureg.day + 4 * ureg.hour + 43 * ureg.minute,
    # 04/16 06:00
    13 * ureg.day + 19 * ureg.hour + 40 * ureg.minute,
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


mass_transport_coeff_factor = 3

k_top = 4.9e-7 * ureg.m * ureg.s**-1 * mass_transport_coeff_factor * 0.16
optimised_ratio = 0.1
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
neutron_rate = (P383_neutron_rate + A325_neutron_rate) * 0.5

baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=calculated_TBR,
    k_top=k_top,
    k_wall=k_wall,
    neutron_rate=neutron_rate,
    irradiations=irradiations,
)
