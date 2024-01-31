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


background = 0.29 * ureg.Bq
vial_511 = background_sub(0.254 * ureg.Bq, background)
vial_512 = background_sub(0.289 * ureg.Bq, background)
vial_513 = background_sub(3.885 * ureg.Bq, background)
vial_514 = background_sub(0.344 * ureg.Bq, background)

vial_521 = background_sub(0.262 * ureg.Bq, background)
vial_522 = background_sub(0.283 * ureg.Bq, background)
vial_523 = background_sub(3.998 * ureg.Bq, background)
vial_524 = background_sub(0.408 * ureg.Bq, background)

vial_531 = background_sub(0.266 * ureg.Bq, background)
vial_532 = background_sub(0.293 * ureg.Bq, background)
vial_533 = background_sub(4.644 * ureg.Bq, background)
vial_534 = background_sub(0.396 * ureg.Bq, background)

vial_541 = background_sub(0.275 * ureg.Bq, background)
vial_542 = background_sub(0.278 * ureg.Bq, background)
vial_543 = background_sub(5.098 * ureg.Bq, background)
vial_544 = background_sub(0.505 * ureg.Bq, background)

vial_551 = background_sub(0.276 * ureg.Bq, background)
vial_552 = background_sub(0.286 * ureg.Bq, background)
vial_553 = background_sub(2.823 * ureg.Bq, background)
vial_554 = background_sub(0.468 * ureg.Bq, background)

vial_561 = background_sub(0.265 * ureg.Bq, background)
vial_562 = background_sub(0.281 * ureg.Bq, background)
vial_563 = background_sub(1.252 * ureg.Bq, background)
vial_564 = background_sub(0.333 * ureg.Bq, background)

vial_571 = background_sub(0.262 * ureg.Bq, background)
vial_572 = background_sub(0.286 * ureg.Bq, background)
vial_573 = background_sub(0.874 * ureg.Bq, background)
vial_574 = background_sub(0.387 * ureg.Bq, background)


baby_diameter = 1.77 * ureg.inches - 2 * 0.06 * ureg.inches  # from CAD drawings
baby_radius = 0.5 * baby_diameter
baby_volume = 0.125 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section
baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=4.57e-4 * ureg.particle * ureg.neutron**-1,  # stefano 1/22/2024
)


mass_transport_coeff_factor = 3 * 0.7

baby_model.k_top *= mass_transport_coeff_factor
baby_model.k_wall *= mass_transport_coeff_factor

baby_model.number_days = 2 * ureg.days
baby_model.exposure_time = 12 * ureg.hour

baby_model.irradiations = [
    [0 * ureg.hour, 0 + baby_model.exposure_time],
    [24 * ureg.hour, 24 * ureg.hour + baby_model.exposure_time],
]

initial_neutron_rate = (
    (1.2e8 + 3.96e8) * ureg.neutron * ureg.s**-1
)  # initially measured by activation foils

fitting_param = 0.82
baby_model.neutron_rate = fitting_param * initial_neutron_rate
baby_model.dt = 0.05 * ureg.h
