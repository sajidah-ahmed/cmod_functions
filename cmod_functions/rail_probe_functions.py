import MDSplus as mds
import numpy as np

variables_dictionary_rail = {
    "Js": "JS_SLOW",
    "Jg": "JGRND_SLOW",
}


def get_rail_signal(
    shot_number,
    variable_name,
    tile,
    probe_pin_number,
    time_start=-np.inf,
    time_end=np.inf,
):

    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    dataname = f"\EDGE::TOP.PROBES.RAIL.TILE{tile}.P{probe_pin_number}:{variables_dictionary_rail[variable_name]}"

    rail_signal = c.get(dataname).data()
    rail_signal_time = c.get(f"dim_of({dataname})").data()

    time_interval = (rail_signal_time > time_start) & (rail_signal_time < time_end)

    return rail_signal_time[time_interval], rail_signal[time_interval]


def get_rail_rho(
    shot_number,
    tile,
    probe_pin_number,
):

    c = mds.Connection("alcdata")
    c.openTree("edge", shot_number)

    dataname_rho = f"\EDGE::TOP.PROBES.RAIL.TILE{tile}.P{probe_pin_number}:RHO"

    rho = c.get(dataname_rho)
    rho_time = c.get(f"dim_of({dataname_rho})").data()

    return rho_time, rho
