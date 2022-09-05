import MDSplus as mds
import numpy as np

"""
Abbreviations:
    FSP: F-port scanning probe (scans in the vertical direction)

    Note: The mirror-Langmuir probe was NOT operational before 2012.

    p0 is the NE pin
    p1 is the SE pin
    p2 is the SW pin
    p3 is the NW pin
"""

# Node names for FSP data. Use this convention.
# Fast data is highly sampled data. Slow data is downsampled data.
variables_dictionary_fsp = {
    "ne_fast": "NE_FAST",
    "Te_fast": "TE_FAST",
    "Is_fast": "I_FAST",
    "Vf_fast": "V_FAST",
    "ne_slow": "NE_SLOW",
    "Te_slow": "TE_SLOW",
    "Is_slow": "I_SLOW",
    "Vf_slow": "V_SLOW",
}


def get_fsp_plunge_depth(shot_number: int):
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

    dataname_plunge = "\EDGE::TOP.PROBES.FSP_1.PLUNGE"

    plunge = c.get(dataname_plunge).data()
    plunge_time = c.get(f"dim_of({dataname_plunge})").data()

    return plunge_time, plunge


def get_fsp_probe_origin(shot_number: int):
    """
    Extracts the probe origin.

    Args:
        shot_number: Shot number(s) of interest.

    Returns:
        origin: Probe origin giving (R, Z, R*Phi) of the probe, *EXPLAIN R, Z, Phi*
    """

    dataname_origin = "\EDGE::TOP.PROBES.FSP_1.G_1.ORIGIN"
    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    origin = c.get(dataname_origin).data()

    return origin


def get_fsp_rho(shot: int, probe_pin_number: int):
    """
    Extracts the rho, the distance relative to the last-closed flux surface, of a FSP probe tip.
    Shots before 2012 have FSP data (i.e. conventional Langmuir probe).

    Args:
        shot_number: Shot number(s) of interest.
        probe_pin_number: Particluar probe tip usually from 0 to 3.

    Returns:
        rho_time: Time data for rho
        rho: The probe position relative to the separatrix
    """

    c = mds.Connection("alcdata")
    c.openTree("edge", shot)

    dataname_rho = f"\EDGE::TOP.PROBES.FSP_1.G_1.P{probe_pin_number}:RHO"

    rho = c.get(dataname_rho)
    rho_time = c.get(f"dim_of({dataname_rho})").data()

    return rho_time, rho


def get_raw_fsp_data(
    shot_number: int,
    probe_pin_number: int,
    variable_name: str,
    time_start: float = -np.inf,
    time_end: float = np.inf,
):
    """
    Extracts raw FSP data. Shots before 2012 have conventional Langmuir probe data.
    This is just a conventional Langmuir probe.

    Args:
        shot_number: Shot number(s) of interest.
        probe_pin_number: Particluar probe tip usually from 0 to 3.
        variable_name: The variable of interests
            variables_dictionary_fsp = {
                "ne_fast": "NE_FAST",
                "Te_fast": "TE_FAST",
                "Is_fast": "I_FAST",
                "Vf_fast": "V_FAST",
                "ne_slow": "NE_SLOW",
                "Te_slow": "TE_SLOW",
                "Is_slow": "I_SLOW",
                "Vf_slow": "V_SLOW",
            }

        time_start: Start time of interest. Set to first frame by default.
        time_end: End time of interest. Set to last frame by default.

    Returns:
        fsp_time: Time data for FSP.
        fsp_data: Raw FSP data of a particular variable.
    """

    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    dataname = f"\EDGE::TOP.PROBES.FSP_1.G_1.P{probe_pin_number}:{variables_dictionary_fsp[variable_name]}"

    fsp_data = c.get(dataname).data()

    fsp_time = c.get(f"dim_of({dataname})").data()

    time_interval = (fsp_time > time_start) & (fsp_time < time_end)
    return fsp_time[time_interval], fsp_data[time_interval]
