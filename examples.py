from cmod_functions import get_major_radius_phantom_coordinates, get_phantom_frames
import xarray as xr

R, Z = get_major_radius_phantom_coordinates(1110114011)

time, frames = get_phantom_frames(1110114011)

ds = xr.Dataset(
        {"frames": (["time", "y", "x"], frames)},
        coords={"R": (["y", "x"], R), "Z": (["y", "x"], Z), "time": (["time"], time)},
    )

ds.to_netcdf('phantom_xarray.nc')
print(ds)