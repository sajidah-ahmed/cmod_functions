import MDSplus as mds
import numpy as np


def get_line_integrated_density(shot_number):
    """
    Extract line integrated density data.
    To get line-averaged density, you can use the function in this package, but it will only exist for certain shots.
    Otherwise divide by the chord length - this is derived from the equilibrium reconstruction.

    Args:
        shot_number: shot number(s) of interest.

    Returns:
        line_integrated_density_time: Time data for the ne_bar data.
        line_integrated_density: Line-integrated density measured in per sqaure metres (m^-2).
    """

    c = mds.Connection("alcdata")
    c.openTree("electrons", shot_number)

    line_integrated_density_dataname = "\ELECTRONS::TOP.TCI.RESULTS.NL_04"

    line_integrated_density = c.get(line_integrated_density_dataname).data()
    line_integrated_density_time = c.get(
        f"dim_of({line_integrated_density_dataname})"
    ).data()

    return line_integrated_density_time, line_integrated_density


def get_line_averaged_density(shot_number):
    """
    Extract line averaged density data.

    Args:
        shot_number: shot number(s) of interest.
                     NOTE: The data for a particular shot may not exist.
                           They're only generated if someone asked for them previously.

    Returns:
        line_averaged_density_time: Time data for the ne_bar data.
        line_averaged_density: Line-averaged density measured in per cubic metres (m^-3).
    """

    c = mds.Connection("alcdata")
    c.openTree("electrons", shot_number)

    line_averaged_density_dataname = "\ELECTRONS::TOP.TCI.RESULTS.INVERSION.NEBAR_EFIT"

    line_averaged_density = c.get(line_averaged_density_dataname).data()
    line_averaged_density_time = c.get(
        f"dim_of({line_averaged_density_dataname})"
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

    plasma_current_dataname = "\MAGNETICS::IP/1000"
    plasma_current = c.get(plasma_current_dataname).data()
    plasma_current_time = c.get(f"dim_of({plasma_current_dataname})").data()

    return plasma_current_time, plasma_current


def get_toroidal_magnetic_field(shot_number):
    """
    Extract toroidal magnetic field data.

    Args:
        shot_number: shot number(s) of interest.

    Returns:
        toroidal_magnetic_field_time: time data for toroidal_magnetic_field data.
        toroidal_magnetic_field: toroidal magnetic field measured in Tesla (T).
                                (-) sign indicates the normal field direction.
                                (+) sign indicates the reverse field direction.
    """

    c = mds.Connection("alcdata")
    c.openTree("magnetics", shot_number)

    toroidal_magnetic_field_dataname = "\MAGNETICS::BTOR"

    toroidal_magnetic_field = c.get(toroidal_magnetic_field_dataname).data()
    toroidal_magnetic_field_time = c.get(
        f"dim_of({toroidal_magnetic_field_dataname})"
    ).data()

    return toroidal_magnetic_field_time, toroidal_magnetic_field

def get_q95(shot_number):
    """
    Extract the safety factor, q at 95% flux surface. This is called q95.

    Args:
        shot_number: shot number(s) of interest.

    Returns:
        q95_time: time data for q95.
        q95: q95 data. This is dimensionless.
    """

    c = mds.Connection("alcdata")
    c.openTree("analysis", shot_number)

    q95_dataname = "\ANALYSIS::EFIT_AEQDSK:QPSIB"

    q95 = c.get(q95_dataname).data()
    q95_time = c.get(
        f"dim_of({q95_dataname})"
    ).data()

    return q95_time, q95


def get_kappa(shot_number):
    """
    Extract the elongation/ellipticity at plasma boundary. This is called kappa.
    This is how it's calculated: http://fusionwiki.ciemat.es/wiki/Ellipticity

    Args:
        shot_number: shot number(s) of interest.

    Returns:
        kappa_time: time data for kappa.
        kappa: kappa data. This is dimensionless.
    """

    c = mds.Connection("alcdata")
    c.openTree("analysis", shot_number)

    kappa_dataname = "\ANALYSIS::EFIT_AEQDSK:EOUT"

    kappa = c.get(kappa_dataname).data()
    kappa_time = c.get(
        f"dim_of({kappa_dataname})"
    ).data()

    return kappa_time, kappa


def get_delta(shot_number,type):
    """
    Extract the traingularity at plasma boundary

    Args:
        shot_number: shot number(s) of interest.
        type: "upper" means you return the upper triangularity
              "lower" means you return the lower triangularity
              "average" means you return the average of the upper and lower triangularity

    Returns:
        delta_time: time data for delta.
        delta: delta data. This is dimensionless.
    """

    c = mds.Connection("alcdata")
    c.openTree("analysis", shot_number)

    if type=="upper":

        delta_dataname = "\ANALYSIS::EFIT_AEQDSK:DOUTU"
        delta = c.get(delta_dataname).data()
        delta_time = c.get(
            f"dim_of({delta_dataname})"
        ).data()
    elif type == "lower":
        
        delta_dataname = "\ANALYSIS::EFIT_AEQDSK:DOUTL"
        delta = c.get(delta_dataname).data()
        delta_time = c.get(
            f"dim_of({delta_dataname})"
        ).data()
    elif type == "average":

        delta_upper_dataname = "\ANALYSIS::EFIT_AEQDSK:DOUTU"
        delta_lower_dataname = "\ANALYSIS::EFIT_AEQDSK:DOUTL"

        delta_upper = c.get(delta_upper_dataname).data()
        delta_lower = c.get(delta_lower_dataname).data()

        delta = (delta_upper + delta_lower) / 2
        delta_time = c.get(
            f"dim_of({delta_upper_dataname})"
        ).data()

    return delta_time, delta


def get_plasma_area(shot_number):
    """
    Extract the poloidal cross-sectional area of the plasma

    Args:
        shot_number: shot number(s) of interest.

    Returns:
        area_time: time data for the area..
        area: area data. This is in units of m^2.
    """

    c = mds.Connection("alcdata")
    c.openTree("analysis", shot_number)

    area_dataname = "\ANALYSIS::EFIT_AEQDSK:AREA"

    area = c.get(area_dataname).data()
    area_time = c.get(
        f"dim_of({area_dataname})"
    ).data()

    return area_time, area


def get_plasma_volume(shot_number):
    """
    Extract the plasma volume

    Args:
        shot_number: shot number(s) of interest.

    Returns:
        volume_time: time data for the volume..
        volume: volume data. This is in units of m^3.
    """

    c = mds.Connection("alcdata")
    c.openTree("analysis", shot_number)

    volume_dataname = "\ANALYSIS::EFIT_AEQDSK:VOLUME"

    volume = c.get(volume_dataname).data()
    volume_time = c.get(
        f"dim_of({volume_dataname})"
    ).data()

    return volume_time, volume


def get_plasma_stored_energy(shot_number):
    """
    Extract the plasma stored energy, in units of Joules (J).

    Args:
        shot_number: shot number(s) of interest.

    Returns:
        wplasma_time: time data for the stored energy.
        wplasma: wplasma data. This is in units of J.
    """

    c = mds.Connection("alcdata")
    c.openTree("analysis", shot_number)

    wplasma_dataname = "\ANALYSIS::EFIT_AEQDSK:WPLASM"

    wplasma = c.get(wplasma_dataname).data()
    wplasma_time = c.get(
        f"dim_of({wplasma_dataname})"
    ).data()

    return wplasma_time, wplasma


def average_plasma_parameter(
    variable_data, variable_time, time_start: float = -np.inf, time_end: float = np.inf
):
    """
    Use: Calculates the average values of your plasma parameters

    Inputs:
        variable_data: Data of your variable.
            This can be the plasma current, the line-averaged density, toroidal magnetic field etc.
        variable_time: Corresponding time data of your variable.
        time_start: The starting time of your time window.
            Default set to the minimum time of time data.
        time_end: The end time of your time window.
            Default set to the maximum time of time data.

    Outputs:
        variable_mean: Mean value of variable in the time window of choice.

    """

    time_interval = (variable_time < time_end) & (variable_time > time_start)
    variable_range = variable_data[time_interval]
    variable_mean = variable_range.mean()

    return variable_mean


def greenwald_density_limit(average_plasma_current, minor_radius=0.22):
    """
    Use: Calculates the Greenwald density limit

    Inputs:
        average_plasma_current: The average plasma current in Mega Amps.
        minor_radius: The minor radius in metres.
                      This is 0.22 metres for Alcator C-Mod by default.

    Outputs:
        greenwald_density_limit: The density limit in units of 10^20 m^-3.

    """
    return average_plasma_current / (np.pi * minor_radius * minor_radius)


def greenwald_fraction(average_line_averaged_density, greenwald_density):
    """
    Use: Calculates the Greenwald fraction which the ratio between the average line-averaged density and the Greenwald density limit.

    Inputs:
        average_line_averaged_density: The average line-averaged density in 10^20 m^-3 units.
        greenwald_density: The Greenwald density limit in 10^20 m^-3 units.

    Outputs:
        greenwald_fraction: The Greenwald fraction - should be between 0 and 1.

    """
    return average_line_averaged_density / greenwald_density
