# %%
import os
import glob
import numpy as np
import scipy.io as io
import utils
import TKED as tked
import gsw
import seawater
from oceans.sw_extras import gamma_GP_from_SP_pt
import utm


data_in = "../raw_data"

glob1 = glob.glob(os.path.join(data_in, "sp14-towyo-*"))
glob2 = glob.glob(os.path.join(data_in, "rr1209_gridded_towyo*"))
towyo_files = np.sort(np.hstack((glob1, glob2)))

TOWYOs = [utils.loadmat(file)["tm"] for file in towyo_files]

# %%
# Sill lat boundaries (very roughly)

P5_lat = -8.5  # Greater than this
P2_lat = -9.5  # Less than this
# P4 is between these

# Reprocessing parameters
acc = 5e-4
R0 = 0.25
save_dir = "../proc_data"
Nsq_method = "bulk"
c1 = 0.9  # The coefficient that multiplies the Thorpe scales... e.g. 0.64, 0.9
sig4b = 1045.9  # The potential density used in the buoyancy estimate 
sig4max = 1045.93
dc = 50  # step size for estimating distances
serr = 34.65  # Lowest salinity that we do not classify as an error

iP5 = 0
iP2 = 0
iP4 = 0
for j, (file, TY) in enumerate(zip(towyo_files, TOWYOs)):
    
    if "sp14" in file:
        TY["year"] = 2014
    else:
        TY["year"] = 2012
    
    # Classify by location and year.
    if np.nanmean(TY["lat"]) > P5_lat:
        TY["sill"] = "P5"
        TY["name"] = "P5_{:02d}".format(iP5)
        TY["longname"] = "TY_" + TY["name"] + "_{}".format(TY["year"])
        iP5 += 1
    elif np.nanmean(TY["lat"]) < P2_lat:
        TY["sill"] = "P2"
        TY["name"] = "P2_{:02d}".format(iP2)
        TY["longname"] = "TY_" + TY["name"] + "_{}".format(TY["year"])
        iP2 += 1
    else:
        TY["sill"] = "P4"
        TY["name"] = "P4_{:02d}".format(iP4)
        TY["longname"] = "TY_" + TY["name"] + "_{}".format(TY["year"])
        iP4 += 1
        
    print("PROCESSING {}".format(TY["longname"]))

    if "eps" in TY:
        TY["eps_gunnar"] = TY.pop("eps")

    if "np" not in TY:
        TY["np"] = len(TY["mlon"])

    _, Np = TY["t1"].shape
    
    print("- Fixing salinity errors by linear inerpolation")
    allbad = TY["s1"] < serr  # will catch nans
    # Don't want to interpolate nans
    nans = np.isnan(TY["s1"])
    errs = allbad.copy()
    errs[nans] = False
    for i in range(Np):
        err = errs[:, i]
        if not err.any():
            continue
            
        good = ~allbad[:, i]
        TY["s1"][err, i] = np.interp(TY["z"][err], TY["z"][good], TY["s1"][good, i])

    print("- Renaming variables and calculating new ones")
    TY["T"] = TY.pop("t1")
    TY["P"] = TY.pop("p")
    TY["S"] = TY.pop("s1")
    TY["spd"] = np.sqrt(TY["u"] ** 2 + TY["v"] ** 2)
    TY["z"] = -np.tile(TY["z"][:, np.newaxis].astype(float), (1, Np))
    TY["depth"] = -TY["z"]
    TY["SA"] = gsw.SA_from_SP(TY["S"], TY["P"], TY["lon"], TY["lat"])
    TY["CT"] = gsw.CT_from_t(TY["SA"], TY["T"], TY["P"])
    TY["sig4"] = gsw.pot_rho_t_exact(TY["SA"], TY["T"], TY["P"], 4000)
    TY["g"] = gsw.grav(TY["lat"], TY["P"])
    TY["b"] = -TY["g"] * (TY["sig4"] - sig4b) / sig4b
    TY["sig4m"] = np.nanmean(TY["sig4"])
    TY["sig4_sorted"] = utils.nansort(TY["sig4"], axis=0)
    TY["b_sorted"] = -TY["g"] * (TY["sig4_sorted"] - sig4b) / sig4b
    TY["N2"], TY["P_mid"] = gsw.Nsquared(
        TY["SA"], TY["CT"], TY["P"], np.nanmean(TY["lat"], axis=0)
    )
    TY["z_mid"] = gsw.z_from_p(TY["P_mid"], TY["mlat"])
    ddist = utils.lldist(TY["mlon"], TY["mlat"])
    TY["dist"] = np.hstack((0, np.cumsum(ddist)))
    TY["dist_r"] = -(TY["dist"] - TY["dist"].max())

    print("- Estimating distances")
    lon = TY["lon"].flatten()
    lat = TY["lat"].flatten()
    time = TY["datenum"].flatten()
    use = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(time)
    idxs = np.argsort(time[use])
    time_ = time[use][idxs]
    lon_ = lon[use][idxs]
    lat_ = lat[use][idxs]
    d = utils.lldist(lon_[::dc], lat_[::dc])
    dint = np.hstack((0, np.cumsum(d)))
    tint = time_[::dc]
    dist_ = np.full_like(lon, np.nan)
    dist_[use] = np.interp(time[use], tint, dint)
    x, y, TY["zone_number"], TY["zone_letter"] = utm.from_latlon(lat_[::dc], lon_[::dc])
    x_, y_ = np.full_like(lon, np.nan), np.full_like(lon, np.nan)
    x_[use] = np.interp(time[use], tint, x)
    y_[use] = np.interp(time[use], tint, y)
    TY["distdata"] = np.reshape(dist_, TY["lon"].shape)
    TY["xdata"] = np.reshape(x_, TY["lon"].shape)
    TY["ydata"] = np.reshape(y_, TY["lon"].shape)
    TY["xm"] = np.nanmean(TY["xdata"], axis=0)
    TY["ym"] = np.nanmean(TY["ydata"], axis=0)
    
    print("- Estimating bottom depth")
    nans = np.isnan(TY["P"])
    depth = TY["depth"].copy()
    depth[nans] = np.nan
    bdepth = np.nanmax(depth, axis=0)
    try:
        bdist = TY["bdist"].copy().astype(float)
        nans = np.isnan(bdist)
        bdist[nans] = 40
    except KeyError:
        bdist = np.full_like(TY["mlon"], 40.0)
        
    TY["bdepth"] = bdepth + bdist

    print("- Calculating the overflow integrated quantities")
    u_ = TY["u"].copy()
    v_ = TY["v"].copy()
    # Note this mask also excludes deep array padding.
    mask = TY["sig4"] < sig4max
    u_[mask] = np.nan
    v_[mask] = np.nan
    TY["uo"] = utils.nantrapz(u_, x=TY["z"], axis=0, xave=True)
    TY["vo"] = utils.nantrapz(v_, x=TY["z"], axis=0, xave=True)
    # Need the minus here because otherwise we get the wrong way around...
    TY["UT"] = -utils.nantrapz(u_, x=TY["z"], axis=0, xave=False)
    TY["VT"] = -utils.nantrapz(v_, x=TY["z"], axis=0, xave=False)
    TY["zo"] = utils.nan_interp(sig4max, TY["sig4_sorted"], TY["z"], axis=0)
    # Density difference 200 m above the interface. 
    TY["dsig4"] = utils.nan_interp(-TY["zo"] - 100., TY["depth"], TY["sig4_sorted"], axis=0)[0, :] - sig4max
    TY["HKE"] = utils.nantrapz(0.5*(u_**2 + v_**2), x=TY["z"], axis=0, xave=True)

    print("- Calculating neutral density")
    TY["PT0"] = seawater.ptmp(TY["S"], TY["T"], TY["P"])
    # Flatten variables for analysis.
    nans = np.isnan(TY["lon"].flatten())
    S = TY["S"].flatten()[~nans]
    T = TY["PT0"].flatten()[~nans]
    P = TY["P"].flatten()[~nans]
    LO = TY["lon"].flatten()[~nans]
    LA = TY["lat"].flatten()[~nans]
    gamman = np.full_like(TY["lon"].flatten(), np.nan)
    gamman[~nans] = gamma_GP_from_SP_pt(S, T, P, LO, LA)
    TY["gamman"] = np.reshape(gamman, TY["S"].shape) + 1000.0
    TY["gamman_sorted"] = utils.nansort(TY["gamman"], axis=0)
    TY["gammanm"] = np.nanmean(TY["gamman"])

    TY["N2_sorted"] = (
        -TY["g"]
        * np.gradient(TY["sig4_sorted"], axis=0)
        / (TY["sig4m"] * np.gradient(TY["z"], axis=0))
    )

    print("- Calculating dissipation rate")
    LT = np.full_like(TY["T"], np.nan)
    N2_overturn = np.full_like(TY["T"], 0.0)
    Lo = np.full_like(TY["T"], 0.0)
    LT_counter = []
    Lo_counter = []
    TY["N2_smoothed"] = np.full_like(TY["N2_sorted"], np.nan)

    for i in range(Np):
        nans = np.isnan(TY["sig4"][:, i])
        if nans.all():
            continue

        z = TY["z"][~nans, i]
        den = TY["sig4"][~nans, i]

        LT_, _, Nsqo, Lo_, _, _, _ = tked.thorpe_scales1(
            z, den, acc=acc, R0=R0, Nsq_method=Nsq_method, full_output=True
        )

        N2_overturn[~nans, i] = Nsqo
        Lo[~nans, i] = Lo_
        LT[~nans, i] = LT_
        idxs = utils.contiguous_regions(LT_ > 0.0)

        for idx in idxs:
            LT_counter.append(LT_[idx[0] : idx[1]][0])
            Lo_counter.append(Lo_[idx[0] : idx[1]][0])

    TY["eps"] = c1 * LT ** 2 * N2_overturn ** 1.5
    TY["N2_overturn"] = N2_overturn
    TY["LT"] = LT
    TY["Lo"] = Lo
    TY["LT_counter"] = LT_counter
    TY["Lo_counter"] = Lo_counter

    nans = np.isnan(TY["eps"])
    eps_ = np.zeros_like(TY["eps"])
    eps_[~nans] = TY["eps"][~nans]
    # I guess this works because we're not averaging...
    TY["eps_zint"] = -1000 * np.trapz(eps_, TY["z"][:, 0], axis=0)
    # This is the density mask... add NaNs back again.
    eps_[mask] = np.nan
    TY["epso"] = utils.nantrapz(eps_, TY["z"][:, 0], axis=0, xave=True)

    io.savemat(os.path.join(save_dir, TY["longname"] + ".mat"), TY)
