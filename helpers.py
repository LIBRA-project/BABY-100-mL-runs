from simple_tritium_transport_model import ureg


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


def substract_background_from_measurements(raw_measurements):
    """Substracts the background from a set of measurements.

    Args:
        raw_measurements (dict): The raw measurements. The keys are the measurement
        numbers, and the values are dictionaries with the keys being the vial
        numbers and the values being the measured activities. The dictionary should
        also contain a key "background" with the background activity.

    Returns:
        dict: The measurements with the background substracted. The keys are the
        measurement numbers, and the values are dictionaries with the keys being the
        vial numbers and the values being the measured activities with the
        background substracted.
    """
    measurements_after_background_sub = {
        i: {
            j: background_sub(act, raw_measurements[i]["background"])
            for j, act in raw_measurements[i].items()
            if j != "background"
        }
        for i in raw_measurements
    }
    return measurements_after_background_sub


def cumulative_activity(measurements):
    """Calculates the cumulative activity of a set of measurements.

    Args:
        measurements (dict): The measurements. The keys are the measurement numbers,
        and the values are dictionaries with the keys being the vial numbers and the
        values being the measured activities.

    Returns:
        list: The cumulative activity of the measurements
    """
    for measurement in measurements.values():
        if "background" in measurement:
            raise ValueError("Background should be substracted first")

    total_samples = []
    for measurement in measurements.values():
        total_sample = sum(list(measurement.values()))
        total_samples.append(total_sample)

    cumulative_values = []
    for total_sample in total_samples:
        if cumulative_values:  # if list is not empty
            cumulative_values.append(cumulative_values[-1] + total_sample)
        else:
            cumulative_values.append(total_sample)

    return cumulative_values
