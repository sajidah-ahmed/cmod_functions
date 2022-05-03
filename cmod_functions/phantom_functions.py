import MDSplus as mds
import numpy as np
import xarray as xr
import scipy as sp
from typing import Union


def get_integration_time(shot_number: int):
    """
    Extract integration time .
    Args:
        shot_number: Shot number of interest.
    Returns:
        integration time
    """
    c = mds.Connection("alcdata")
    c.openTree("spectroscopy", shot_number)

    return c.get("GPI.PHANTOM.SETTINGS:EXPOS").data()


def get_sensitivity_callibration(shot_number: int, filter: str):
    """
    Extract sensitivity callibration data.
    Args:
        shot_number: Shot number of interest.
        filter: H_alpha or He-587 filter
    Returns:
        inverse sensitivity for every pixel
    """
    assert filter in {"H_alpha", "He-587"}

    c = mds.Connection("alcdata")
    c.openTree("spectroscopy", shot_number)

    if filter == "H_alpha":
        return c.get("GPI.PHANTOM:SENS_ARR_DA").data()
    if filter == "He-587":
        return c.get("GPI.PHANTOM:SENS_ARR_HEI").data()


def get_limiter_coordinates(shot_number: int):
    """
    Extracts the radial and poloidal positions of the limiter shadow in major radius coordinates for shots.
    Args:
        shot_number: Shot number of interest.
    Returns:
        R_limiter: Major radius coordinates in centimetres.
        Z_limiter: Height array (above the machine midplane) in centimetres.
    """
    c = mds.Connection("alcdata")
    c.openTree("mhd", shot_number)
    R_limiter = c.get("\MHD::TOP.analysis.limiters.gh_limiter:R")
    Z_limiter = c.get("\MHD::TOP.analysis.limiters.gh_limiter:Z")

    return R_limiter, Z_limiter


def get_separatrix_coordinates(shot_number: int):
    """
    Extracts the radial and poloidal positions of the last closed flux surface (LCFS) in major radius coordinates for shots.
    Args:
        shot_number: Shot number(s) of interest.
    Returns:
        R_LCFS: Major radius coordinates in centimetres.
        Z_LCFS: Height array (above the machine midplane) in centimetres.
        time_LCFS: time array

    """
    c = mds.Connection("alcdata")
    c.openTree("analysis", shot_number)
    rbbbs = c.get("\efit_geqdsk:rbbbs")
    zbbbs = c.get("\efit_geqdsk:zbbbs")
    nbbbs = c.get("\efit_geqdsk:nbbbs")
    efit_time = c.get("dim_of(\efit_geqdsk:rbbbs)")

    return rbbbs, zbbbs, nbbbs, efit_time


def calculate_splinted_LCFS(
    time_step: float,
    efit_time: np.ndarray,
    rbbbs: np.ndarray,
    zbbbs: np.ndarray,
    nbbbs: np.ndarray,
):
    time_index = np.where(
        np.abs(efit_time - time_step) == np.min(np.abs(efit_time - time_step))
    )
    z_fine = np.linspace(-0.06, 0.01, 100)
    r_index = np.where(rbbbs[: nbbbs[time_index[0]], time_index[0]] >= 0.86)
    r_sep = rbbbs[r_index[:][0], time_index[0]]
    z_sep = zbbbs[r_index[:][0], time_index[0]]

    zsep_sort = np.sort(z_sep)
    jt = np.array(np.argsort(z_sep))
    rsep_sort = r_sep[jt]

    f = sp.interpolate.CubicSpline(zsep_sort, rsep_sort, bc_type="natural")
    r_fine = f(z_fine)

    return r_fine, z_fine


def get_major_radius_phantom_coordinates(shot_number: int):
    """
    Extracts the radial and poloidal positions of the pixels in major radius coordinates for shots.
    Args:
        shot_number: Shot number(s) of interest.
    Returns:
        R: Major radius coordinates in centimetres.
        Z: Height array (above the machine midplane) in centimetres.
    """

    c = mds.Connection("alcdata")
    c.openTree("spectroscopy", shot_number)

    R_array_path = "GPI.PHANTOM.IMAGE_POS.R_ARR"
    Z_array_path = "GPI.PHANTOM.IMAGE_POS.Z_ARR"

    R = c.get(R_array_path).data()
    Z = c.get(Z_array_path).data()

    return R, Z


def get_phantom_frames(shot_number: int):
    """
    Extracts the frames and time array for a shot.

    Args:
        shot_number: Shot number(s) of interest.

    Returns:
        times: Time array.
        frames: Raw frames extracted for all pixels.
    """

    c = mds.Connection("alcdata")
    c.openTree("spectroscopy", shot_number)
    frames_path = "GPI.PHANTOM.FRAMES"
    frame_rate_path = "GPI.PHANTOM.SETTINGS.FRAME_RATE"
    t_hists_path = "GPI.PHANTOM.SETTINGS.TRIG_TIME"

    frames = c.get(frames_path).data()

    frame_rate = c.get(frame_rate_path).data()
    start_time = c.get(t_hists_path).data()

    number_of_frames = frames.shape[0]
    stop_time = start_time + number_of_frames / frame_rate

    times = np.arange(start=start_time, stop=stop_time, step=1 / frame_rate)

    return times, frames


def _normalize_frames(
    frames: np.ndarray,
    shot_number: int,
    filter: Union[None, int] = None,
):
    """
    !!!     Check normalization with Jim    !!!
    """
    try:
        inverse_sensitivity = get_sensitivity_callibration(
            shot_number=shot_number, filter=filter
        )
        integration_time = get_integration_time(shot_number)
        return np.matmul(frames, inverse_sensitivity) * integration_time
    except Exception:
        print(f"No sensitivity callibration data available")
        return frames


def generate_phantom_dataset(
    shot_number: int,
    time_start: float = -np.inf,
    time_end: float = np.inf,
    filter: Union[None, str] = None,
):
    """
    Generates an xarray dataset containing raw outboard midplane Phantom camera data for a shot

    Args:
        shot_number: Shot number(s) of interest.
        time_start: The beginning of the time window in seconds. Set to first frame by default.
        time_end: The end of the time window in seconds. Set to last frame by default.
        filter: if not None, using callibration filter

    Returns:
        dataset: An xarray dataset containing raw phantom camera data for all pixels:
            x, y = [0,0] refers to the lowest R and Z values.

            time: Time array in seconds.
            frames: Raw frames extracted for all pixels.
            R: Major radius coordinates in centimetres if available.
            Z: Height array (above the machine midplane) in centimetres if avaliable.
    """
    assert filter in {None, "H_alpha", "He-587"}

    time, frames = get_phantom_frames(shot_number)
    frames = np.flip(frames, axis=2)

    time_interval = (time > time_start) & (time < time_end)
    frames = frames[time_interval, :, :]
    time = time[time_interval]

    if filter is not None:
        frames = _normalize_frames(frames, shot_number, filter)

    # Not all shots provide R and Z data
    try:
        R, Z = get_major_radius_phantom_coordinates(shot_number)
        R = np.flip(R, axis=1)
        Z = np.flip(Z, axis=1)

        return xr.Dataset(
            {"frames": (["time", "y", "x"], frames)},
            coords={
                "R": (["y", "x"], R),
                "Z": (["y", "x"], Z),
                "time": (["time"], time),
            },
        )

    except Exception:
        print(f"No R and Z data available for shot {shot_number}")

        return xr.Dataset(
            {"frames": (["time", "y", "x"], frames)},
            coords={"time": (["time"], time)},
        )
