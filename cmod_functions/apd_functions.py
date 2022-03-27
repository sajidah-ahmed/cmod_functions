import MDSplus as mds
import numpy as np
import shot_details


def get_major_radius_coordinates(shot_number):
    """
    Extracts the radial and poloidal positions of the pixels in major radius coordinates for shots.

    Args:
        shot_number: Shot number(s) of interest.

    Returns:
        R_array: Radial coordinates in centimetres.
        Z_array: Poloidal coordinates in centrimetres.
    """

    c = mds.Connection("alcdata")
    c.openTree("spectroscopy", shot_number)

    R_array_path = "GPI.APD_ARRAY.R_ARR"
    Z_array_path = "GPI.APD_ARRAY.Z_ARR"

    # extract frames, time, radial array and poloidal array
    R_array = c.get(R_array_path).data()
    Z_array = c.get(Z_array_path).data()

    return R_array, Z_array


def get_raw_apd_frames(shot_number):
    """
    Extracts the frames and time array for a shot.

    Args:
        shot_number: Shot number(s) of interest.

    Returns:
        times_array: Time array in 100 nanosecond units.
        frames_array: Frames extracted for all pixels.
    """

    c = mds.Connection("alcdata")
    c.openTree("spectroscopy", shot_number)
    frames_path = "GPI.APD_ARRAY.FRAMES"

    frames_array = c.get(frames_path).data()
    time_array = c.get("dim_of(" + frames_path + ")").data()

    return time_array, frames_array


def generate_apd_data(shot_number, location):
    """
    Generate the frames and time array of all the pixels for a shot. You must make a shot_details.py file with
    the time window of interest for each shot.

    Args:
        shot_number: Shot number(s) of interest.
        location: Locations of pixel(s).

    Returns:
        times_array: Time array in seconds.
        signal_array: Unnormalized time series extracted for pixels in the time window of interest.
    """

    time_array, frames_array = get_raw_apd_frames(shot_number)

    # Invert voltage
    frames = -frames_array

    # Convert to seconds
    time = time_array[-frames.shape[0] :] * 1e-7

    # Sometimes, the time dimension is not the same as the number of frames.
    # need to force them to be the same by adding time points if the time dimension is smaller than the number of frames
    # or cutting the time values off the end of the time array if the time dimension is greater than the number of frames.

    time_array_length = len(time_array)
    frame_array_length = len(frames_array)

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

    time_start, time_end = shot_details.apd_time_dictionary[shot_number]
    time_interval = (time > time_start) & (time < time_end)
    frames = frames[time_interval, :, :]
    time = time[time_interval]

    # Choose frames
    location_array_length = len(location)
    moving_window = 8192

    # Create new array of smaller size by chopping one window size on both ends
    signal_array = np.zeros(
        (
            location_array_length,
            frames[2 * moving_window : -2 * moving_window, 0, 0].size,
        )
    )

    time_array = time[2 * moving_window : -2 * moving_window]

    for i in range(len(location)):
        raw_time_series = frames[:, location[i][0], location[i][1]]
        signal_array[i, :] = raw_time_series[:]

    # Reshape signal array
    signal_array = np.swapaxes(np.reshape(signal_array, (10, 9, len(time_array))), 0, 1)

    return time_array, signal_array


def efit_major_radius_to_rho(R, Z, time_array, shot_number, tree):
    """
    Converts radial and poloidal coordinates to flux surface coordinates, rho.

    Args:
        R: Major radius (in cm) array for the coordinate pair to be mapped to rho.
        Z: Height (above the machine midplane) (in cm) array for the coordinate pair to be mapped to rho.
        time_array: Time array in seconds for the times for which the mapping is desired.
        shot_number: Shot number of interest.
        tree: Set this equal to the string of the tree name (e.g. 'efit19') which will be used
                for the flux surface mapping. The default is 'analysis'

    Returns:
        rho: 2D array (with the 1st dimension having the same number of
            elements as inputs R and Z and the 2nd dimension having the
            same number of elements as the input time_array).
            The 1st dimension specifies the element of the coordinate pair (R, Z) for the
            flux surface mapping from the coordinate (R, Z) at time time_array(i)
            to the distance beyond the separatrix in the horizontal plane
            that is defined by the Z of the magnetic axis (ZMAGX) at
            time time_array(i). Rho is returned in cm, with positive sign being outside
            the separatrix and minus sign being inside the separatrix.

    """

    import eqtools as eq
    import numpy as np

    cmod_efit = eq.CModEFITTree(shot_number, tree)

    # Major radii (in m) of the R,Z coordinate arrays
    # flux-surface-mapped to the height of the magnetic axis at each
    # time of the time_array. NOTE: that the time index is the 1st one, the space index is the 2nd!
    rmid = cmod_efit.rz2rmid(R / 100.0, Z / 100.0, time_array)

    # rmidout is the major radius of the outboard side of the separatrix at the
    # height of the magnetic axis (m)
    rmidout = cmod_efit.getRmidOut()
    efit_time = cmod_efit.getTimeBase()

    rho = np.zeros((len(R), len(time_array)))

    for i, ti in enumerate(time_array[:]):

        # Rhos is now the array of distances between the R,Z coordinate array
        # points that have been flux-surface-mapped to the height of the
        # magnetic axis MINUS the major radius of the outboard side of the
        # sepx at the height of the magnetic axis (in cm)

        conditional = np.where(
            np.abs(time_array[i] - efit_time)
            == np.amin(np.abs(time_array[i] - efit_time))
        )
        rho[:, i] = (rmid[i, :] - rmidout[conditional]) * 100.0

    return rho


def major_radius_to_rho(shot_number, location):

    R_array, Z_array = get_major_radius_coordinates(shot_number)

    time_start, time_end = shot_details.apd_time_dictionary[shot_number]

    # Map R, Z major radius coordinates onto a flux surface
    # Rho is the distance from the last-closed flux surface in centimetres
    location_array_length = len(location)
    time = np.arange(time_start, time_end, 0.001)

    # Make a list of the R, Z coordinates
    major_radius_position = []

    for i in range(location_array_length):

        R = R_array[location[i][0], location[i][1]]
        Z = Z_array[location[i][0], location[i][1]]

        RZ_pair = [R, Z]
        major_radius_position.append(RZ_pair)

    major_radius_R = np.array(major_radius_position)[:, 0]
    major_radius_Z = np.array(major_radius_position)[:, 1]

    rho_array = efit_major_radius_to_rho(
        R=major_radius_R,
        Z=major_radius_Z,
        time_array=time,
        shot_number=shot_number,
        tree="EFIT19",
    )

    rho_mean = np.mean(rho_array, axis=1)
    time_averaged_rho = np.swapaxes(np.reshape(rho_mean, (10, 9)), 0, 1)

    # Reshape array
    major_radius_R_array = np.swapaxes(np.reshape(major_radius_R, (10, 9)), 0, 1)
    major_radius_Z_array = np.swapaxes(np.reshape(major_radius_Z, (10, 9)), 0, 1)

    return major_radius_R_array, major_radius_Z_array, time_averaged_rho
