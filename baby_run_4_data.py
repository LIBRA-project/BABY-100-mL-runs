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
