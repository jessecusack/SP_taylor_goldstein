import numpy as np
import glob
import munch
import utils
import xarray as xr


def stack1Dvar(TYs, var):
    return np.hstack([TY[var] for TY in TYs])


def stack2Dvar(TYs, var, zmin=None, zmax=None, dz=1., npfl=None):
    if zmin is None:
        zmin = np.min([TY.z.min() for TY in TYs])
    if zmax is None:
        zmax = np.max([TY.z.max() for TY in TYs])
    if npfls is None:
        npfl = np.sum([TY.mlon.shape[0] for TY in TYs_])
        
    z = np.arange(zmax, zmin + dz, dz)
    nz = len(z)
    
    # Initialise stacked array
    vara = np.full((nz, npfl), np.nan)
    j0 = 0
    for TY in TYs:
        z0 = TY.z[0, 0]
        z1 = TY.z[-1, 0]
        j1 = j0 + TY.mlon.shape[0]  # Profile index
        i0, i1 = np.searchsorted(-z, [-z0, -z1])  # Minus needed because z array is monotonically decreasing.
        
        # Fill stacked array, the +1 is needed because of the way searchsorted returns indices on the left.
        vara[i0:i1+1, j0:j1] = TY[var]
        j0 = j1
        
    return vara


if __name__ == "__main__":
    # Load processed towyos.
    data_dir = "../proc_data"
    data_files = glob.glob(data_dir + "/TY*.mat")
    TYs = np.asarray([munch.munchify(utils.loadmat(file)) for file in data_files])

    # Get P5 towyos from 2014.
    TYs_ = [TY for TY in TYs if ((TY.sill == "P5") and (TY.year == 2014))]

    zmin = np.min([TY.z.min() for TY in TYs_])
    zmax = np.max([TY.z.max() for TY in TYs_])
    npfls = np.sum([TY.mlon.shape[0] for TY in TYs_])

    idxa = np.arange(npfls)
    dz = -1.
    z = np.arange(zmax, zmin + dz, dz)
    nz = len(z)

    datavars = {}

    # Fill data dictionary
    varlist = ["T", "S", "CT", "SA", "u", "v", "sig4", "sig4_sorted", "b", "eps", "LT", "Lo", "N2_overturn"]
    for var in varlist:
        datavars[var] = (["z", "pfl"], stack2Dvar(TYs_, var, zmin, zmax, dz, npfls))


    # Fill 1D coords
    coords = {
        "z": (["z"], z),
        "pfl": (["pfl"], idxa),
        "lon": (["pfl"], stack1Dvar(TYs_, "mlon")),
        "lat": (["pfl"], stack1Dvar(TYs_, "mlat")),
        "x": (["pfl"], stack1Dvar(TYs_, "xm")),
        "y": (["pfl"], stack1Dvar(TYs_, "ym")),
    }

    varlist = ["mlon", "mlat"]

    ds = xr.Dataset(datavars, coords)

    # Save to netcdf.
    ds.to_netcdf(data_dir + "/stacked_towyos.nc")