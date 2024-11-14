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

file_reader_3 = LSCFileReader(
    "data/IV_1-3-X_BL-1.csv",
    vial_labels=[
        "IV-BL-1",
        None,
        "IV-1-3-1",
        "IV-1-3-2",  # probably has a statistic issue
        "IV-1-3-3",
        "IV-1-3-4",
    ],
)
file_reader_3.read_file()

file_reader_4 = LSCFileReader(
    "data/Report1.csv",
    vial_labels=[
        "BL-1_count_4",
        None,
        "IV-1-4-1",
        "IV-1-4-2",
        "IV-1-4-3",
        "IV-1-4-4",
        None,
        "OV-1-1-1",
        "OV-1-1-2",
        "OV-1-1-3",
        "OV-1-1-4",
        None,
        "IV-1-3-2 (repeat)",
    ],
)
file_reader_4.read_file()


# Make samples

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

sample_3_IV = LIBRASample(
    samples=[
        LSCSample.from_file(file_reader_3, "IV-1-3-1"),
        LSCSample.from_file(
            file_reader_4, "IV-1-3-2 (repeat)"
        ),  # the first one has a statistic issue
        LSCSample.from_file(file_reader_3, "IV-1-3-3"),
        LSCSample.from_file(file_reader_3, "IV-1-3-4"),
    ],
    time="11/10/2024 1:33 PM",
)
blank_sample_3_IV = LSCSample.from_file(file_reader_3, "IV-BL-1")


sample_4_IV = LIBRASample(
    samples=[
        LSCSample.from_file(file_reader_4, label)
        for label in ["IV-1-4-1", "IV-1-4-2", "IV-1-4-3", "IV-1-4-4"]
    ],
    time="11/13/2024 2:31 PM",
)
blank_sample_4 = LSCSample.from_file(file_reader_4, "BL-1_count_4")

sample_1_OV = LIBRASample(
    samples=[
        LSCSample.from_file(file_reader_4, label)
        for label in ["OV-1-1-1", "OV-1-1-2", "OV-1-1-3", "OV-1-1-4"]
    ],
    time="11/13/2024 2:31 PM",
)

# Make streams

start_time = "11/4/2024 10:07 AM"

IV_stream = GasStream(
    [sample_1_IV, sample_2_IV, sample_3_IV, sample_4_IV], start_time=start_time
)
OV_stream = GasStream([sample_1_OV], start_time=start_time)

# substract background
for sample in [sample_1_IV, sample_2_IV]:
    sample.substract_background(
        background_sample=LSCSample(activity=0.320 * ureg.Bq, name="background")
    )  # TODO don't have a real background here

sample_3_IV.substract_background(background_sample=blank_sample_3_IV)
sample_4_IV.substract_background(background_sample=blank_sample_4)
sample_1_OV.substract_background(background_sample=blank_sample_4)

# create run
run = LIBRARun(streams=[IV_stream, OV_stream], start_time=start_time)

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
replacement_times_walls = [
    sample.get_relative_time(start_time) for sample in OV_stream.samples
]

# convert timedelta to pint quantity  # TODO add this to libra-toolbox

replacement_times_top = [
    time.total_seconds() * ureg.second for time in replacement_times_top
]
replacement_times_walls = [
    time.total_seconds() * ureg.second for time in replacement_times_walls
]

replacement_times_top = sorted(replacement_times_top)

replacement_times_walls = sorted(replacement_times_walls)


# Model

baby_diameter = 14 * ureg.cm  # TODO confirm with CAD
baby_radius = 0.5 * baby_diameter
baby_volume = 1 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section

# from OpenMC
calculated_TBR = 2.0e-3 * ureg.particle * ureg.neutron**-1

optimised_ratio = 3e-2
k_top = 7.8e-8 * ureg.m * ureg.s**-1
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
