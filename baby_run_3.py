from simple_tritium_transport_model import ureg, Model
import numpy as np
from helpers import substract_background_from_measurements

background_1 = 0.265 * ureg.Bq

raw_measurements = {
    1: {
        1: 0.279 * ureg.Bq,
        2: 0.284 * ureg.Bq,
        3: 6.352 * ureg.Bq,
        4: 0.617 * ureg.Bq,
        "background": background_1,
    },
    2: {
        1: 0.272 * ureg.Bq,
        2: 0.278 * ureg.Bq,
        3: 6.621 * ureg.Bq,
        4: 0.581 * ureg.Bq,
        "background": background_1,
    },
    3: {
        1: 0.252 * ureg.Bq,
        2: 0.272 * ureg.Bq,
        3: 2.347 * ureg.Bq,
        4: 0.706 * ureg.Bq,
        "background": background_1,
    },
}

measurements_after_background_sub = substract_background_from_measurements(
    raw_measurements
)

replacement_times = [
    1 * ureg.day,
    2 * ureg.day,
    4 * ureg.day,
]

m = measurements_after_background_sub
sample_1 = sum(list(m[1].values()))
sample_2 = sum(list(m[2].values()))
sample_3 = sum(list(m[3].values()))

cumulative_1 = sample_1
cumulative_2 = sample_1 + sample_2
cumulative_3 = sample_1 + sample_2 + sample_3

cumulative_values = [cumulative_1, cumulative_2, cumulative_3]

baby_diameter = 1.77 * ureg.inches - 2 * 0.06 * ureg.inches  # from CAD drawings
baby_radius = 0.5 * baby_diameter
baby_volume = 0.125 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section
baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=3.3e-4 * ureg.particle * ureg.neutron**-1,  # stefano 10/24/2023
)

fitting_param = 0.86

mass_transport_coeff_factor = 3

baby_model.k_top *= mass_transport_coeff_factor
baby_model.k_wall *= mass_transport_coeff_factor

baby_model.number_days = 2 * ureg.days
baby_model.exposure_time = 12 * ureg.hour
baby_model.neutron_rate = fitting_param * (1.2e8 + 3.96e8) * ureg.neutron * ureg.s**-1
baby_model.dt = 0.4 * ureg.h
