import MDSplus as mds
import numpy as np
import xarray as xr

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

    R = c.get(R_array_path)#.data()
    Z = c.get(Z_array_path)#.data()

    return R, Z

def get_phantom_frames(shot_number: int):
    """
    Extracts the frames and time array for a shot.

    Args:
        shot_number: Shot number(s) of interest.

    Returns:
        times: Time array.
        frames: Raw frames extracted for all pixels. These require further processing before analysis.
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

    stop_time = start_time + number_of_frames/frame_rate

    times = np.arange(start = start_time, stop = stop_time, step = 1/frame_rate)

    return times, frames

