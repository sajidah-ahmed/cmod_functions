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


def get_apd_frames(shot_number):
    """
    Extracts the frames and time array for a shot.

    Args:
        shot_number: Shot number(s) of interest.

    Returns:
        times_array: Time array in 100 nanosecond units.
        frames_array: Raw frames extracted for all pixels. These require further processing before analysis.
    """

    c = mds.Connection("alcdata")
    c.openTree("spectroscopy", shot_number)
    frames_path = "GPI.APD_ARRAY.FRAMES"

    frames_array = c.get(frames_path).data()
    time_array = c.get("dim_of(" + frames_path + ")").data()

    return time_array, frames_array


def generate_raw_apd_dataset(shot_number):
    """
    Generates an xarray dataset containing raw APD data for a shot

    Args:
        shot_number: Shot number(s) of interest.

    Returns:
        dataset: An xarray dataset containing raw APD data for all pixels:
            time: Time array in 100 nanosecond units.
            frames: Raw frames extracted for all pixels. These require further processing before analysis.
            R: Major radius coordinates in centimetres.
            Z: Height array (above the machine midplane) in centimetres.
    """

    time, frames = get_apd_frames(shot_number)
    R, Z = get_major_radius_coordinates(shot_number)

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

    import xarray as xr

    frames = np.swapaxes(np.reshape(frames, (9, 10, len(time))), 0, 1)

    dataset = xr.Dataset(
        {"frames": (["y", "x", "time"], frames)},
        coords={"R": (["y", "x"], R), "Z": (["y", "x"], Z), "time": (["time"], time)},
    )

    return dataset


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


def major_radius_to_average_rho(shot_number, location, time_slice=False, tree="EFIT19"):
    """
    Given the pixel locations and the time slice, this function converts radial and poloidal coordinates to flux surface coordinates, rho.

    Args:
        shot_number: Shot number of interest.
        location: Pixel location where pixel = (Z, R). Z is the height (above the machine midplane) and R is the major radius, both in centimetres.
        time_slice: Time window of the signal to be analysed. Default is 'False' which takes the time range of the APD switched on.
                    If set 'True' make a shot_details.py script with time windows specified.
        tree: Set this equal to the string of the tree name (e.g. 'efit19') which will be used
                for the flux surface mapping. The default is 'analysis'

    Returns:
        major_radius_R_array: Array of radial coordinates in centimetres.
        major_radius_Z_array: Array of the height in centimetres.
        time_averaged_rho: The time-averaged flux mapped coordinate in centimetres based on the time window specified.

    """

    R_array, Z_array = get_major_radius_coordinates(shot_number)

    if time_slice:
        import shot_details

        time_start, time_end = shot_details.apd_time_dictionary[shot_number]
    else:
        time_array, _ = get_apd_frames(shot_number)
        time_start, time_end = time_array.amin(), time_array.amax()

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
        tree=tree,
    )

    rho_mean = np.mean(rho_array, axis=1)
    time_averaged_rho = np.swapaxes(np.reshape(rho_mean, (9, 10)), 0, 1)

    # Reshape array
    major_radius_R_array = np.swapaxes(np.reshape(major_radius_R, (9, 10)), 0, 1)
    major_radius_Z_array = np.swapaxes(np.reshape(major_radius_Z, (9, 10)), 0, 1)

    return major_radius_R_array, major_radius_Z_array, time_averaged_rho
