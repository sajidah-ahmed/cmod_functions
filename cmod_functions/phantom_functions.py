import MDSplus as mds
import numpy as np
import xarray as xr


def get_limiter_coordinates(shot_number: int):
    """
    Extracts the radial and poloidal positions of the limiter shadow in major radius coordinates for shots.
    Args:
        shot_number: Shot number(s) of interest.
    Returns:
        R: Major radius coordinates in centimetres.
        Z: Height array (above the machine midplane) in centimetres.
    """
    c = mds.Connection("alcdata")
    c.openTree("mhd", shot_number)
    R_limiter = c.get("\MHD::TOP.analysis.limiters.gh_limiter:R")
    Z_limiter = c.get("\MHD::TOP.analysis.limiters.gh_limiter:Z")

    return R_limiter, Z_limiter


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


def generate_phantom_dataset(
    shot_number: int, time_start: float = -np.inf, time_end: float = np.inf
):
    """
    Generates an xarray dataset containing raw outboard midplane Phantom camera data for a shot

    Args:
        shot_number: Shot number(s) of interest.
        time_start: The beginning of the time window in seconds. Set to first frame by default.
        time_end: The end of the time window in seconds. Set to last frame by default.

    Returns:
        dataset: An xarray dataset containing raw phantom camera data for all pixels:
            x, y = [0,0] refers to the lowest R and Z values.

            time: Time array in seconds.
            frames: Raw frames extracted for all pixels.
            R: Major radius coordinates in centimetres if available.
            Z: Height array (above the machine midplane) in centimetres if avaliable.
    """

    time, frames = get_phantom_frames(shot_number)
    frames = np.flip(frames, axis=2)

    time_interval = (time > time_start) & (time < time_end)
    frames = frames[time_interval, :, :]
    time = time[time_interval]

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
