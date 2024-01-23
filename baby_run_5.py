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


# time starts at 12/05 9:30 AM
# 12/06 9:30 AM = 24 hours
# 12/07 9:30 AM = 48 hours
replacement_times = [
    # 12/05 22:58
    0 * ureg.day + 13 * ureg.hour + 28 * ureg.minute,
    # 12/06 09:09
    0 * ureg.day + 23 * ureg.hour + 39 * ureg.minute,
    # 12/06 21:58
    1 * ureg.day + 12 * ureg.hour + 28 * ureg.minute,
    # 12/07 11:19
    1 * ureg.day + 25 * ureg.hour + 49 * ureg.minute,
    # 12/08 09:53
    2 * ureg.day + 24 * ureg.hour + 23 * ureg.minute,
    # 12/09 08:28
    3 * ureg.day + 22 * ureg.hour + 58 * ureg.minute,
    # 12/11 12:48
    6 * ureg.day + 3 * ureg.hour + 18 * ureg.minute,
]

replacement_times = sorted(replacement_times)

baby_diameter = 1.77 * ureg.inches - 2 * 0.06 * ureg.inches  # from CAD drawings
baby_radius = 0.5 * baby_diameter
baby_volume = 0.085 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section
baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=4.78e-4 * ureg.particle * ureg.neutron**-1,  # stefano 1/22/2024
)


mass_transport_coeff_factor = 3 * 0.6 * 0.9

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

fitting_param = 0.72
baby_model.neutron_rate = fitting_param * initial_neutron_rate
baby_model.dt = 0.05 * ureg.h
