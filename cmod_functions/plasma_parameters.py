import MDSplus as mds


def get_line_averaged_density(shot_number):
    """
    Extract line averaged density data.

    Args:
        shot_number: shot number(s) of interest.

    Returns:
        line_averaged_density_time: Time data for the ne_bar data.
        line_averaged_density: Line-averaged density measured in per cubic metres (m^-3).
    """

    c = mds.Connection("alcdata")
    c.openTree("electrons", shot_number)

    line_averaged_density_dataname = f"\ELECTRONS::TOP.TCI.RESULTS.INVERSION.NEBAR_EFIT"

    line_averaged_density = c.get(line_averaged_density_dataname).data()
    line_averaged_density_time = c.get(
        "dim_of(" + line_averaged_density_dataname + ")"
    ).data()

    return line_averaged_density_time, line_averaged_density


def get_plasma_current(shot_number):
    """
    Extract plasma current trace

    Args:
        shot_number: shot number(s) of interest.

    Returns:
        plasma_current_time: Time data for the plasma current data.
        plasma_current: Plasma current measures in  kilo Amps.
                        Negative sign means it's in the normal field direction.
                        Positive sign means it's in the reverse field direction.
    """

    c = mds.Connection("alcdata")
    c.openTree("magnetics", shot_number)

    plasma_current_dataname = f"\MAGNETICS::IP/1000"
    plasma_current = c.get(plasma_current_dataname).data()
    plasma_current_time = c.get("dim_of(" + plasma_current_dataname + ")").data()

    return plasma_current_time, plasma_current


def get_toroidal_magnetic_field(shot_number):
    """
    Extract toroidal magnetic field data.

    Args:
        shot_number: shot number(s) of interest.

    Returns:
        toroidal_magnetic_field_time: time data for toroidal_magnetic_field data.
        toroidal_magnetic_field: toroidal magnetic field measured in Tesla (T).
    """

    c = mds.Connection("alcdata")
    c.openTree("magnetics", shot_number)

    toroidal_magnetic_field_dataname = f"\MAGNETICS::BTOR"

    toroidal_magnetic_field = -c.get(toroidal_magnetic_field_dataname).data()
    toroidal_magnetic_field_time = c.get(
        "dim_of(" + toroidal_magnetic_field_dataname + ")"
    ).data()

    return toroidal_magnetic_field_time, toroidal_magnetic_field
