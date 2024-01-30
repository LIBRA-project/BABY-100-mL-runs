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


background = 0.2308 * ureg.Bq
vial_611 = background_sub(0.2388 * ureg.Bq, background)
vial_612 = background_sub(0.2442 * ureg.Bq, background)
vial_613 = background_sub(3.5140 * ureg.Bq, background)
vial_614 = background_sub(0.2950 * ureg.Bq, background)

vial_621 = background_sub(0.2599 * ureg.Bq, background)
vial_622 = background_sub(0.3967 * ureg.Bq, background)
vial_623 = background_sub(2.4481 * ureg.Bq, background)
vial_624 = background_sub(0.7118 * ureg.Bq, background)

# time starts at 01/25 9:36 AM
# 01/26 9:36 AM = 24 hours
# 01/27 9:36 AM = 48 hours
replacement_times = [
    # 01/25 21:44
    0 * ureg.day + 15 * ureg.hour + 8 * ureg.minute,
    # 01/26 09:22
    1 * ureg.day + 0 * ureg.hour - 14 * ureg.minute,
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
