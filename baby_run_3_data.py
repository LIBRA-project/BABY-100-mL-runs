from simple_tritium_transport_model import ureg, Model
import numpy as np


def background_sub(measured, background):
    """Substracts the background of a measured activity.
    Returns zero if the background is greater than measurement.

    Args:
        measured (pint.Quantity): The measured activity
        background (pint.Quantity): the background acitivity

    Returns:
        pint.Quantity: activity with substracted background
    """
    if measured > background:
        return measured - background
    else:
        return 0 * ureg.Bq


# day 1 sample
background_1 = 0.265 * ureg.Bq
vial_311 = background_sub(0.279 * ureg.Bq, background_1)
vial_312 = background_sub(0.284 * ureg.Bq, background_1)
vial_313 = background_sub(6.352 * ureg.Bq, background_1)
vial_314 = background_sub(0.617 * ureg.Bq, background_1)


# day 2 sample
background_2 = 0.274 * ureg.Bq
vial_321 = background_sub(0.272 * ureg.Bq, background_2)
vial_322 = background_sub(0.278 * ureg.Bq, background_2)
vial_323 = background_sub(6.621 * ureg.Bq, background_2)
vial_324 = background_sub(0.581 * ureg.Bq, background_2)

# day 4 sample
background_3 = 0.252 * ureg.Bq
vial_331 = background_sub(0.252 * ureg.Bq, background_3)
vial_332 = background_sub(0.272 * ureg.Bq, background_3)
vial_333 = background_sub(2.347 * ureg.Bq, background_3)
vial_334 = background_sub(0.706 * ureg.Bq, background_3)

baby_diameter = 1.77 * ureg.inches - 2 * 0.06 * ureg.inches  # from CAD drawings
baby_radius = 0.5 * baby_diameter
baby_volume = 0.1 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section
baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=3.3e-4 * ureg.particle * ureg.neutron**-1,  # stefano 10/24/2023
)

fitting_param = 0.82

mass_transport_coeff_factor = 3

baby_model.k_top *= mass_transport_coeff_factor
baby_model.k_wall *= mass_transport_coeff_factor

baby_model.number_days = 2 * ureg.days
baby_model.exposure_time = 12 * ureg.hour
baby_model.neutron_rate = fitting_param * (1.2e8 + 3.96e8) * ureg.neutron * ureg.s**-1
baby_model.dt = 0.4 * ureg.h
