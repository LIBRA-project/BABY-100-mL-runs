from libra_toolbox.tritium.model import ureg, Model
import numpy as np
from libra_toolbox.tritium.helpers import (
    substract_background_from_measurements,
    cumulative_activity,
    background_sub,
)


raw_measurements_top = {}

measurements_top_after_background_sub = substract_background_from_measurements(
    raw_measurements_top
)


# time starts at 11/4/2024 10:00 AM
replacement_times_top = [
    # 11/5/2024 10:00 AM
    24
    * ureg.hour,
]

replacement_times_top = sorted(replacement_times_top)

replacement_times_walls = []

replacement_times_walls = sorted(replacement_times_walls)

# # Cumulative values

cumulative_release_top = cumulative_activity(measurements_top_after_background_sub)


# Model

baby_diameter = 14 * ureg.cm  # TODO confirm with CAD
baby_radius = 0.5 * baby_diameter
baby_volume = 1 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section

# from OpenMC
calculated_TBR = 2.5e-3 * ureg.particle * ureg.neutron**-1
baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=calculated_TBR,
)

mass_transport_coeff_factor = 3

baby_model.k_top *= mass_transport_coeff_factor * 0.7
optimised_ratio = 3e-2
baby_model.k_wall = baby_model.k_top * optimised_ratio

exposure_time = 12 * ureg.hour

baby_model.irradiations = [
    [0 * ureg.hour, 0 + exposure_time],
]

# calculated from Kevin's activation foil analysis
P383_neutron_rate = 4.95e8 * ureg.neutron * ureg.s**-1
A325_neutron_rate = 2.13e8 * ureg.neutron * ureg.s**-1

neutron_rate_relative_uncertainty = 0.089
baby_model.neutron_rate = (
    A325_neutron_rate
) / 2  # the neutron rate is divided by two to acount for the double counting (two detectors)
