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


background_2 = 0.28 * ureg.Bq
vial_421 = background_sub(0.299 * ureg.Bq, background_2)
vial_422 = background_sub(0.277 * ureg.Bq, background_2)
vial_423 = background_sub(0.283 * ureg.Bq, background_2)
vial_424 = background_sub(0.274 * ureg.Bq, background_2)

background_3 = background_2  # no sure what the blank is for 3
vial_431 = background_sub(0.279 * ureg.Bq, background_3)
vial_432 = background_sub(0.276 * ureg.Bq, background_3)
vial_433 = background_sub(0.564 * ureg.Bq, background_3)
vial_434 = background_sub(0.290 * ureg.Bq, background_3)

background_4 = background_2  # no sure what the blank is for 4
vial_441 = background_sub(0.331 * ureg.Bq, background_4)
vial_442 = background_sub(0.296 * ureg.Bq, background_4)
vial_443 = background_sub(2.100 * ureg.Bq, background_4)
vial_444 = background_sub(0.310 * ureg.Bq, background_4)

background_5 = background_2  # no sure what the blank is for 5
vial_451 = background_sub(0.269 * ureg.Bq, background_5)
vial_452 = background_sub(0.284 * ureg.Bq, background_5)
vial_453 = background_sub(4.050 * ureg.Bq, background_5)
vial_454 = background_sub(0.369 * ureg.Bq, background_5)


background_6 = background_2  # no sure what the blank is for 6
vial_461 = background_sub(0.247 * ureg.Bq, background_6)
vial_462 = background_sub(0.308 * ureg.Bq, background_6)
vial_463 = background_sub(8.469 * ureg.Bq, background_6)
vial_464 = background_sub(0.754 * ureg.Bq, background_6)

background_7 = background_2  # no sure what the blank is for 7
vial_471 = background_sub(0.280 * ureg.Bq, background_7)
vial_472 = background_sub(0.292 * ureg.Bq, background_7)
vial_473 = background_sub(3.439 * ureg.Bq, background_7)
vial_474 = background_sub(0.553 * ureg.Bq, background_7)

background_8 = background_2  # no sure what the blank is for 8
vial_481 = background_sub(0.260 * ureg.Bq, background_8)
vial_482 = background_sub(0.275 * ureg.Bq, background_8)
vial_483 = background_sub(1.243 * ureg.Bq, background_8)
vial_484 = background_sub(0.480 * ureg.Bq, background_8)

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

fitting_param = 1.1

mass_transport_coeff_factor = 3

baby_model.k_top *= mass_transport_coeff_factor
baby_model.k_wall *= mass_transport_coeff_factor

baby_model.number_days = 2 * ureg.days
baby_model.exposure_time = 12 * ureg.hour

baby_model.irradiations = [
    [0 * ureg.hour, 0 + baby_model.exposure_time],
    [24 * ureg.hour, 24 * ureg.hour + baby_model.exposure_time],
]

baby_model.neutron_rate = fitting_param * (1.2e8 + 3.96e8) * ureg.neutron * ureg.s**-1
baby_model.dt = 0.4 * ureg.h
