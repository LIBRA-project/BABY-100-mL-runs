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


background = 0.314 * ureg.Bq
vial_611 = background_sub(0.288 * ureg.Bq, background)
vial_612 = background_sub(0.310 * ureg.Bq, background)
vial_613 = background_sub(4.455 * ureg.Bq, background)
vial_614 = background_sub(0.396 * ureg.Bq, background)

vial_621 = background_sub(0.306 * ureg.Bq, background)
vial_622 = background_sub(0.300 * ureg.Bq, background)
vial_623 = background_sub(3.127 * ureg.Bq, background)
vial_624 = background_sub(0.409 * ureg.Bq, background)

background = 0.301 * ureg.Bq
vial_631 = background_sub(0.358 * ureg.Bq, background)
vial_632 = background_sub(0.286 * ureg.Bq, background)
vial_633 = background_sub(5.735 * ureg.Bq, background)
vial_634 = background_sub(0.442 * ureg.Bq, background)

vial_641 = background_sub(0.272 * ureg.Bq, background)
vial_642 = background_sub(0.299 * ureg.Bq, background)
vial_643 = background_sub(3.915 * ureg.Bq, background)
vial_644 = background_sub(0.398 * ureg.Bq, background)

vial_651 = background_sub(0.301 * ureg.Bq, background)
vial_652 = background_sub(0.305 * ureg.Bq, background)
vial_653 = background_sub(2.930 * ureg.Bq, background)
vial_654 = background_sub(0 * ureg.Bq, background)  # missing water!!

vial_661 = background_sub(0.301 * ureg.Bq, background)
vial_662 = background_sub(0.495 * ureg.Bq, background)
vial_663 = background_sub(1.191 * ureg.Bq, background)
vial_664 = background_sub(0.361 * ureg.Bq, background)

# time starts at 01/25 9:36 AM
# 01/26 9:36 AM = 24 hours = 1 * ureg.day + 0 * ureg.hour + 0 * ureg.minute
# 01/27 9:36 AM = 48 hours = 2 * ureg.day + 0 * ureg.hour + 0 * ureg.minute
replacement_times = [
    # 01/25 21:44
    0 * ureg.day + 15 * ureg.hour + 8 * ureg.minute,
    # 01/26 09:22
    1 * ureg.day
    + (9 * ureg.hour + 22 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
    # 01/26 21:46
    1 * ureg.day
    + (21 * ureg.hour + 46 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
    # 01/27 10:48
    2 * ureg.day
    + (10 * ureg.hour + 48 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
    # 01/28 11:19
    3 * ureg.day
    + (11 * ureg.hour + 19 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
    # 01/29 08:02
    4 * ureg.day
    + (8 * ureg.hour + 2 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
    # 01/30 17:00
    # 5 * ureg.day
    # + (17 * ureg.hour + 0 * ureg.minute)
    # - (9 * ureg.hour + 36 * ureg.minute),
]

replacement_times = sorted(replacement_times)

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
