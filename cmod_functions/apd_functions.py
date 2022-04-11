import MDSplus as mds
import numpy as np
import xarray as xr


def get_major_radius_coordinates(shot_number: int):
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

    R_array_path = "GPI.APD_ARRAY.R_ARR"
    Z_array_path = "GPI.APD_ARRAY.Z_ARR"

    # extract frames, time, radial array and poloidal array
    R = c.get(R_array_path).data()
    Z = c.get(Z_array_path).data()

    return R, Z


def get_apd_frames(shot_number: int):
    """
    Extracts the frames and time array for a shot.

    Args:
        shot_number: Shot number(s) of interest.

    Returns:
        times: Time array in 100 nanosecond units.
        frames: Raw frames extracted for all pixels. These require further processing before analysis.
    """

    c = mds.Connection("alcdata")
    c.openTree("spectroscopy", shot_number)
    frames_path = "GPI.APD_ARRAY.FRAMES"

    frames = c.get(frames_path).data()
    time = c.get("dim_of(" + frames_path + ")").data()

    # Sometimes, the time dimension is not the same as the number of frames.
    # need to force them to be the same by adding time points if the time dimension is smaller than the number of frames
    # or cutting the time values off the end of the time array if the time dimension is greater than the number of frames.

    time_array_length = len(time)
    frame_array_length = len(frames)

    if time_array_length < frame_array_length:
        time_step = (time[99] - time[0]) / 99.0
        num = frame_array_length - time_array_length + 1
        new_times = np.arange(
            time[time_array_length - 1] + time_step,
            time[time_array_length - 1] + num * time_step,
            time_step,
        )
        time = np.append(time, new_times)
    elif time_array_length > frame_array_length:
        time = time[0:frame_array_length]

    return time, frames


def generate_raw_apd_dataset(
    shot_number: int, time_start=None, time_end=None, subtract_background=False
):
    """
    Generates an xarray dataset containing raw APD data for a shot

    Args:
        shot_number: Shot number(s) of interest.
        time_start: The beginning of the time window in seconds. Set to None by default.
        time_end: The end of the time window in seconds. Set to None by default.
        subtract_background: Option to subtract low light levels which will return an inverted signal.
                            Default to False, where this will return just the raw signal, uninverted.

    Returns:
        dataset: An xarray dataset containing raw APD data for all pixels:
            time: Time array in seconds.
            frames: Raw frames extracted for all pixels. These require further processing before analysis.
            R: Major radius coordinates in centimetres.
            Z: Height array (above the machine midplane) in centimetres.
    """

    time, frames = get_apd_frames(shot_number)

    # Convert time from 100 nanosecond units to seconds
    time = time * 1e-7

    # Each pixel may not have the same time dimension as each other.
    # Chop one window size on both ends
    moving_window = 4196  # Equivalent to 2ms of data
    time = time[2 * moving_window : -2 * moving_window]

    R, Z = get_major_radius_coordinates(shot_number)

    apd_signal_array = _create_apd_signal_array(
        frames, moving_window, subtract_background
    )
    return _create_xr_dataset(apd_signal_array, time, time_start, time_end, R, Z)


def _create_apd_signal_array(frames, moving_window, subtract_background):
    """
    Creates an APD signal array from the raw APD frames.
    Time series of dead pixels are replaces with nan values.

    Args:
        frames: Raw frames extracted for all pixels.
        moving_window: The size of the moving window.
        subtract_background: Option to subtract low light levels which will return an inverted signal.

    Returns:
        apd_signal_array: An APD signal array from the raw APD frames.
    """

    apd_pixel_list = np.zeros((90, 2))
    for i in range(90):
        apd_pixel_list[i] = (i % 10, int(i / 10))

    apd_pixel_list = apd_pixel_list.astype(int).tolist()
    pixel_array_length = len(apd_pixel_list)

    apd_signal_array = np.zeros(
        (pixel_array_length, frames[2 * moving_window : -2 * moving_window, 0, 0].size)
    )

    for i in range(len(apd_pixel_list)):
        raw_signal = frames[:, apd_pixel_list[i][0], apd_pixel_list[i][1]]

        # Criterion to find dead pixels
        if raw_signal.std() < 0.01:
            raw_signal[:] = np.nan
        else:
            if subtract_background:
                offset = np.mean(raw_signal[:200])
                raw_signal = offset - raw_signal[:]

        time_series = raw_signal[2 * moving_window : -2 * moving_window]
        apd_signal_array[i, :] = time_series[:]

    return apd_signal_array


def _create_xr_dataset(apd_signal_array, time, time_start, time_end, R, Z):
    """
    Creates an xarray dataset from the raw APD signal array.

    Args:
        apd_signal_array: Raw APD signal array.
        time: Time array in seconds.
        time_start: The beginning of the time window in seconds. Set to None by default.
        time_end: The end of the time window in seconds. Set to None by default.
        R: Major radius coordinates in centimetres.
        Z: Height array (above the machine midplane) in centimetres.

    Returns:
        dataset: An xarray dataset containing raw APD data for all pixels:
            time: Time array in seconds.
            frames: Raw frames extracted for all pixels. These require further processing before analysis.
            R: Major radius coordinates in centimetres.
            Z: Height array (above the machine midplane) in centimetres.
    """

    frames = np.swapaxes(np.reshape(apd_signal_array, (9, 10, len(time))), 0, 1)

    if (time_start is not None) & (time_end is not None):
        time_interval = (time > time_start) & (time < time_end)
        frames = frames[:, :, time_interval]
        time = time[time_interval]

    return xr.Dataset(
        {"frames": (["y", "x", "time"], frames)},
        coords={"R": (["y", "x"], R), "Z": (["y", "x"], Z), "time": (["time"], time)},
    )
