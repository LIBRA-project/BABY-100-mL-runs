from libra_toolbox.tritium.model import ureg, Model
import numpy as np
from libra_toolbox.tritium.lsc_measurements import (
    LSCFileReader,
    LIBRASample,
    LIBRARun,
    LSCSample,
    GasStream,
)

# read files
file_reader_1 = LSCFileReader(
    "data/1L_OV_IV_1-0-X_IV_1-1-X.csv",
    vial_labels=[
        "OV-1-0-1",
        "OV-1-0-2",
        "OV-1-0-3",
        "OV-1-0-4",
        None,
        "IV-1-0-1",
        "IV-1-0-2",
        "IV-1-0-3",
        "IV-1-0-4",
        None,
        "IV-1-1-1",
        "IV-1-1-2",
        "IV-1-1-3",
        "IV-1-1-4",
    ],
)
file_reader_1.read_file()

file_reader_2 = LSCFileReader(
    "data/1L_IV_1-2-X.csv",
    vial_labels=[
        "IV-1-2-1",
        "IV-1-2-2",
        "IV-1-2-3",
        "IV-1-2-4",
    ],
)
file_reader_2.read_file()


sample_0_IV = LIBRASample(
    samples=[
        LSCSample.from_file(file_reader_1, label)
        for label in ["IV-1-0-1", "IV-1-0-2", "IV-1-0-3", "IV-1-0-4"]
    ],
    time="before run",
)

sample_1_IV = LIBRASample(
    samples=[
        LSCSample.from_file(file_reader_1, label)
        for label in ["IV-1-1-1", "IV-1-1-2", "IV-1-1-3", "IV-1-1-4"]
    ],
    time="11/5/2024 12:32 PM",
)

sample_2_IV = LIBRASample(
    samples=[
        LSCSample.from_file(file_reader_2, label)
        for label in ["IV-1-2-1", "IV-1-2-2", "IV-1-2-3", "IV-1-2-4"]
    ],
    time="11/7/2024 8:49 AM",
)

start_time = "11/4/2024 10:07 AM"

IV_stream = GasStream([sample_1_IV, sample_2_IV], start_time=start_time)

# substract background
for sample in IV_stream.samples:
    sample.substract_background(
        background_sample=LSCSample(activity=0.320 * ureg.Bq, name="background")
    )  # TODO don't have a real background here


# create run
run = LIBRARun(streams=[IV_stream], start_time=start_time)

# check that background is always substracted
for stream in run.streams:
    for sample in stream.samples:
        for lsc_vial in sample.samples:
            assert (
                lsc_vial.background_substracted
            ), f"Background not substracted for {sample}"


replacement_times_top = [
    sample.get_relative_time(start_time) for sample in IV_stream.samples
]

# convert timedelta to pint quantity  # TODO add this to libra-toolbox

replacement_times_top = [
    time.total_seconds() * ureg.second for time in replacement_times_top
]

replacement_times_top = sorted(replacement_times_top)

replacement_times_walls = []

replacement_times_walls = sorted(replacement_times_walls)


# Model

baby_diameter = 14 * ureg.cm  # TODO confirm with CAD
baby_radius = 0.5 * baby_diameter
baby_volume = 1 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section

# from OpenMC
calculated_TBR = 2.5e-3 * ureg.particle * ureg.neutron**-1

optimised_ratio = 3e-2
k_top = 7.2e-8 * ureg.m * ureg.s**-1
k_wall = optimised_ratio * k_top

exposure_time = 12 * ureg.hour

irradiations = [
    [0 * ureg.hour, 0 + exposure_time],
]

# calculated from Kevin's activation foil analysis from run 100 mL #7
# TODO replace for values for this run
P383_neutron_rate = 4.95e8 * ureg.neutron * ureg.s**-1
A325_neutron_rate = 2.13e8 * ureg.neutron * ureg.s**-1

neutron_rate_relative_uncertainty = 0.089
neutron_rate = (
    A325_neutron_rate
) / 2  # the neutron rate is divided by two to acount for the double counting (two detectors)

baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=calculated_TBR,
    neutron_rate=neutron_rate,
    irradiations=irradiations,
    k_top=k_top,
    k_wall=k_wall,
)
