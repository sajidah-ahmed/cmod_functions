import MDSplus as mds
import numpy as np

"""
Abbreviations:
    MLP: mirror-Langmuir probe
    ASP: A-port scanning probe (scans in the horizonal direction)
    ISP: Ion-sensitive probe
    
    Functions saying: asp_mlp or asp_isp means they are NOT a convential Langmuir probe
    Note: The mirror-Langmuir probe was NOT operational before 2012.
"""

# Node names for ASP MLP data. Use this convention.
variables_dictionary_asp_mlp = {
    "ne": "DENSITY_FIT",
    "Is": "ISAT_FIT",
    "Js": "JSAT_FIT",
    "Vp": "PHI_FIT",
    "Te": "TE_FIT",
    "Vf": "VF_FIT",
}

# Node names for ASP FAST data. Use this convention.
variables_dictionary_asp = {
    "Is_fast": "I_FAST",
    "Vf_fast": "V_FAST",
    "Is_slow": "I_SLOW",
    "Vf_slow": "V_SLOW",
}

# Node names for ASP ISP data. Use this convention.
variables_dictionary_asp_isp = {
    "Is_fast": "I_FAST",
    "Vf_fast": "V_FAST",
    "Is_slow": "I_SLOW",
    "Vf_slow": "V_SLOW",
}


def get_plunge_depth(shot_number: int):
    """
    Extracts the plunge distance and corresponding time of the probe.

    Args:
        shot_number: Shot number(s) of interest.

    Returns:
        plunge_time: Corresponding time data in seconds.
        plunge: Plunge depth of probe in metres and assumed purely horizontal.
    """

    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    dataname_plunge = "\EDGE::TOP.PROBES.ASP.PLUNGE"

    plunge = c.get(dataname_plunge).data()
    plunge_time = c.get(f"dim_of({dataname_plunge})").data()

    return plunge_time, plunge


def get_probe_origin(shot_number: int):
    """
    Extracts the probe origin.

    Args:
        shot_number: Shot number(s) of interest.

    Returns:
        origin: Probe origin giving (R, Z, R*Phi) of the probe, *EXPLAIN R, Z, Phi*
    """

    dataname_origin = "\EDGE::TOP.PROBES.ASP.G_1.ORIGIN"
    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    origin = c.get(dataname_origin).data()

    return origin


def get_asp_rho(shot: int, probe_pin_number: int):
    """
    Extracts the rho, the distance relative to the last-closed flux surface, of a ASP probe tip.
    Shots before 2012 have ASP data (i.e. conventional Langmuir probe).

    Args:
        shot_number: Shot number(s) of interest.
        probe_pin_number: Particluar probe tip usually from 0 to 3.

    Returns:
        rho_time: Time data for rho
        rho: The probe position relative to the separatrix
    """

    c = mds.Connection("alcdata")
    c.openTree("edge", shot)

    dataname_rho = f"\EDGE::TOP.PROBES.ASP.P{probe_pin_number}:RHO"

    rho = c.get(dataname_rho)
    rho_time = c.get(f"dim_of({dataname_rho})").data()

    return rho_time, rho


def get_raw_asp_data(
    shot_number: int,
    probe_pin_number: int,
    variable_name: str,
    time_start: float = -np.inf,
    time_end: float = np.inf,
):
    """
    Extracts raw ASP data. Shots before 2012 have ASP data.
    This is just a conventional Langmuir probe.

    Args:
        shot_number: Shot number(s) of interest.
        probe_pin_number: Particluar probe tip usually from 0 to 3.
        variable_name: The variable of interests
            variables_dictionary_asp = {
            "Is_fast": "I_FAST",
            "Vf_fast": "V_FAST",
            "Is_slow": "I_SLOW",
            "Vf_slow": "V_SLOW",}
        time_start: Start time of interest. Set to first frame by default.
        time_end: End time of interest. Set to last frame by default.

    Returns:
        asp_time: Time data for ASP.
        asp_data: Raw ASP data of a particular variable.
    """

    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    dataname = f"\EDGE::TOP.PROBES.ASP.G_1.P{probe_pin_number}:{variables_dictionary_asp[variable_name]}"

    asp_data = c.get(dataname).data()

    asp_time = c.get(f"dim_of({dataname})").data()

    time_interval = (asp_time > time_start) & (asp_time < time_end)
    return asp_time[time_interval], asp_data[time_interval]


def get_asp_isp_rho(shot_number: int, probe_pin_number: int):
    """
    Extracts the rho, the distance relative to the last-closed flux surface, of ISP probe tip.
    Check the logbook whether the shot you're after used the ISP.

    Args:
        shot_number: Shot number(s) of interest.
        probe_pin_number: Particluar probe tip usually from 0 to 3.

    Returns:
        rho_time: Time data for rho
        rho: The probe position relative to the separatrix
    """

    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    dataname_rho = f"\EDGE::TOP.PROBES.ASP.ISP.P{probe_pin_number}:RHO"

    rho = c.get(dataname_rho)
    rho_time = c.get(f"dim_of({dataname_rho})").data()

    return rho_time, rho


def get_raw_asp_isp_data(
    shot_number: int,
    probe_pin_number: int,
    variable_name: str,
    time_start: float = -np.inf,
    time_end: float = np.inf,
):
    """
    Extracts raw ISP data. Check the logbook whether the shot you're after used the ISP.

    Args:
        shot_number: Shot number(s) of interest.
        probe_pin_number: Particluar probe tip usually from 0 to 3.
        variable_name: The variable of interests
            variables_dictionary_asp_isp = {
            "Is_fast": "I_FAST",
            "Vf_fast": "V_FAST",
            "Is_slow": "I_SLOW",
            "Vf_slow": "V_SLOW",}
        time_start: Start time of interest. Set to first frame by default.
        time_end: End time of interest. Set to last frame by default.

    Returns:
        asp_isp_time: Time data for ISP.
        asp__ispdata: Raw ISP data of a particular variable.
    """

    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    dataname = f"\EDGE::TOP.PROBES.ASP.ISP.P{probe_pin_number}:{variables_dictionary_asp_isp[variable_name]}"

    asp_isp_data = c.get(dataname).data()

    asp_isp_time = c.get(f"dim_of({dataname})").data()

    time_interval = (asp_isp_time > time_start) & (asp_isp_time < time_end)
    return asp_isp_time[time_interval], asp_isp_data[time_interval]


def get_asp_mlp_rho(shot_number: int, probe_pin_number: int):
    """
    Extracts the rho, the distance relative to the last-closed flux surface, of an MLP probe tip.
    Shots from 2012 onwards have MLP data.

    Args:
        shot_number: Shot number(s) of interest.
        probe_pin_number: Particluar probe tip usually from 0 to 3.

    Returns:
        rho_time: Time data for rho
        rho: The probe position relative to the separatrix
    """

    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    dataname_rho = f"\EDGE::TOP.PROBES.ASP.MLP.P{probe_pin_number}:RHO"

    rho = c.get(dataname_rho)
    rho_time = c.get(f"dim_of({dataname_rho})").data()

    return rho_time, rho


def get_raw_asp_mlp_data(
    shot_number: int,
    probe_pin_number: int,
    variable_name: str,
    time_start: float = -np.inf,
    time_end: float = np.inf,
):
    """
    Extracts raw mirror-Langmuir probe (MLP) data. Shots from 2012 onwards have MLP data. Please interrogate the logbooks.

    Args:
        shot_number: Shot number(s) of interest.
        probe_pin_number: Particluar probe tip usually from 0 to 3.
        variable_name: The variable of interests
            variables_dictionary_asp_mlp = {
            "ne": "DENSITY_FIT",
            "Is": "ISAT_FIT",
            "Js": "JSAT_FIT",
            "Vp": "PHI_FIT",
            "Te": "TE_FIT",
            "Vf": "VF_FIT"}
        time_start: Start time of interest. Set to first frame by default.
        time_end: End time of interest. Set to last frame by default.

    Returns:
        asp_mlp_time: Time data for MLP.
        asp_mlp_data: Raw MLP data of a particular variable.
    """

    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    dataname = f"\EDGE::TOP.PROBES.ASP.MLP.P{probe_pin_number}:{variables_dictionary_asp_mlp[variable_name]}"

    asp_mlp_data = c.get(dataname).data()

    asp_mlp_time = c.get(f"dim_of({dataname})").data()

    time_interval = (asp_mlp_time > time_start) & (asp_mlp_time < time_end)
    return asp_mlp_time[time_interval], asp_mlp_data[time_interval]


def generate_average_mlp_data(shot_number: int, variable_name: str):
    """
    Generates average raw MLP data.

    Args:
        shot_number: Shot number(s) of interest.
        probe_pin_number: Particluar probe tip usually from 0 to 3.
        variable_name: The variable of interests
            variables_dictionary_asp_mlp = {
            "ne": "DENSITY_FIT",
            "Is": "ISAT_FIT",
            "Js": "JSAT_FIT",
            "Vp": "PHI_FIT",
            "Te": "TE_FIT",
            "Vf": "VF_FIT"}

    Returns:
        time_common: Common time data of all four probe bins collecting data.
        time_series_average: Average raw MLP data of a particular variable from all four pin.
    """

    from functools import reduce

    time_series_list = []
    time_list = []

    for probe_pin_number in [0, 1, 2, 3]:
        asp_mlp_time, asp_mlp_data = get_raw_asp_mlp_data(
            variable_name, probe_pin_number, shot_number
        )

        time_series_list.append(asp_mlp_data)
        time_list.append(asp_mlp_time)

    time_common = reduce(np.intersect1d, time_list)

    time_series_average = 0.25 * (
        time_series_list[0]
        + time_series_list[1]
        + time_series_list[2]
        + time_series_list[3]
    )

    return (
        time_common,
        time_series_average,
    )
