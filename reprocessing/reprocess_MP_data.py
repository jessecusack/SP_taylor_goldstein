import os
import glob
import numpy as np
import scipy.io as io
import utils
import TKED as tked
import gsw
import seawater
from oceans.sw_extras import gamma_GP_from_SP_pt

M_dir = '../raw_data/m*_original_data.mat'
Mmoorings_paths = np.sort(glob.glob(M_dir))
P_dir = '../raw_data/p*_original_data.mat'
Pmoorings_paths = np.sort(glob.glob(P_dir))
T_dir = '../raw_data/t*_original_data.mat'
Tmoorings_paths = np.sort(glob.glob(T_dir))

Mmoorings = [utils.loadmat(file) for file in Mmoorings_paths]
Pmoorings = [utils.loadmat(file) for file in Pmoorings_paths]
Tmoorings = [utils.loadmat(file) for file in Tmoorings_paths]
for array in [Mmoorings, Pmoorings, Tmoorings]:
    for i, m in enumerate(array):
        keys = m.keys()
        k = [k for k in keys if '__' not in k][0]
        array[i] = m[k]

MMPs = [Mmoorings[4]['mp'], Mmoorings[5]['mp']]
PMPs = [m['mp'] for m in Pmoorings if m['name'] != 'P2']
TMPs = [m['mp'] for m in Tmoorings]
MPs = MMPs + PMPs + TMPs

altMPs = utils.loadmat('../raw_data/sp_process_moorings_mp.mat',
                       check_arrays=True)['mp']


# Define some convenience functions
def mean_flow_angle(u, v):
    um = np.nanmean(u)
    vm = np.nanmean(v)
    return np.arctan2(vm, um)


def nan_data(mask, x):
    x_ = x.copy()
    x_[mask] = np.nan
    return x_


def rotate_overflow(u, v, rho, rhomax):
    mask = rho < rhomax
    u_ = nan_data(mask, u)
    v_ = nan_data(mask, v)
    a = mean_flow_angle(u_, v_)
    return utils.rotate(u, v, a)


# Reprocessing parameters
acc = 0.5e-3
R0 = 0.25
win = 50
save_dir = '../proc_data'
sig4max = 1045.95
Nsq_method = 'bulk'
c1 = 0.9  # The coefficient that multiplies the Thorpe scales... e.g. 0.64, 0.9

# Do the reprocessing
for MP in MPs:
    if 'name' in MP:
        name = MP['name']
    elif 'info' in MP:
        name = MP['info']['station']
        MP['name'] = name

    print('PROCESSING {}'.format(name))

    if 'eps' in MP:
        MP['eps_gunnar'] = MP.pop('eps')

    # Get gunnar's eps estimate and stick it in here too.
    for aMP in altMPs:
        if aMP['name'] == name:
            MP['eps_gunnar'] = aMP['eps']
            print("- Adding Gunnar's epsilon estimate")

    # Where datenum is NaN we remove those columns from all data.
    nnans = ~np.isnan(MP['datenum'])
    if nnans.all():
        pass
    else:
        print('- Eliminating data columns where datenum is NaN')
        for key in MP:
            item = MP[key]
            if type(item) is not np.ndarray:
                continue
            ndim = np.ndim(item)
            if ndim == 1 and item.size == nnans.size:
                MP[key] = item[nnans]
            elif ndim == 2 and item.shape[1] == nnans.size:
                MP[key] = item[:, nnans]

    print('- Renaming variables and calculating new ones')
    Np = len(MP['id'])
    if np.ndim(MP['lon']) != 0:
        MP['lon'] = MP.pop('lon')[0]
        MP['lat'] = MP.pop('lat')[0]
    MP['P'] = np.tile(MP['p'][:, np.newaxis].astype(float), (1, Np))
    MP['T'] = MP.pop('t')
    MP['S'] = MP.pop('s')
    MP['z'] = gsw.z_from_p(MP['P'], MP['lat'])
    MP['HAB'] = MP['z'] - np.nanmin(MP['z']) + 20.
    MP['SA'] = gsw.SA_from_SP(MP['S'], MP['P'], MP['lon'], MP['lat'])
    MP['CT'] = gsw.CT_from_t(MP['SA'], MP['T'], MP['P'])
    MP['sig4'] = gsw.pot_rho_t_exact(MP['SA'], MP['T'], MP['P'], 4000)
    MP['g'] = gsw.grav(MP['lat'], MP['P'])
    MP['sig4m'] = np.nanmean(MP['sig4'])
    MP['sig4_sorted'] = utils.nansort(MP['sig4'], axis=0)
    MP['N2'], MP['P_mid'] = gsw.Nsquared(MP['SA'], MP['CT'], MP['P'], MP['lat'])
    MP['z_mid'] = gsw.z_from_p(MP['P_mid'], MP['lat'])
    MP['ual'], MP['uac'] = rotate_overflow(MP['u'], MP['v'], MP['sig4'], sig4max)

    print('- Calculating the overflow integrated velocity')
    ual_ = MP['ual'].copy()
    uac_ = MP['uac'].copy()
    mask = MP['sig4'] < sig4max
    ual_[mask] = np.nan
    uac_[mask] = np.nan
    MP['uoal'] = utils.nantrapz(ual_, x=MP['z'], axis=0, xave=True)
    MP['uoac'] = utils.nantrapz(uac_, x=MP['z'], axis=0, xave=True)

    print('- Calculating neutral density')
    MP['PT0'] = seawater.ptmp(MP['S'], MP['T'], MP['P'])
    # Flatten variables for analysis.
    S = MP['S'].flatten()
    T = MP['PT0'].flatten()
    P = MP['P'].flatten()
    LO = MP['lon']*np.ones_like(S)
    LA = MP['lat']*np.ones_like(S)
    gamman = gamma_GP_from_SP_pt(S, T, P, LO, LA)
    MP['gamman'] = np.reshape(gamman, MP['S'].shape) + 1000.
    MP['gamman_sorted'] = utils.nansort(MP['gamman'], axis=0)
    MP['gammanm'] = np.nanmean(MP['gamman'])

    MP['N2_sorted'] = -MP['g']*np.gradient(MP['sig4_sorted'], axis=0) \
        / (MP['sig4m']*np.gradient(MP['z'], axis=0))

    print('- Calculating dissipation rate')
    LT = np.full_like(MP['T'], np.nan)
    N2_overturn = np.full_like(MP['T'], 0.)
    Lo = np.full_like(MP['T'], 0.)
    LT_counter = []
    Lo_counter = []
    MP['N2_smoothed'] = np.full_like(MP['N2_sorted'], np.nan)

    for i in range(Np):
        nans = np.isnan(MP['sig4'][:, i])
        if nans.all():
            continue

        z = MP['z'][~nans, i]
        den = MP['sig4'][~nans, i]

        LT_, _, Nsqo, Lo_, _, _, _ = \
            tked.thorpe_scales1(z, den, acc=acc, R0=R0, Nsq_method=Nsq_method,
                                full_output=True)

        N2_overturn[~nans, i] = Nsqo
        Lo[~nans, i] = Lo_
        LT[~nans, i] = LT_
        idxs = utils.contiguous_regions(LT_ > 0.)

        for idx in idxs:
            LT_counter.append(LT_[idx[0]:idx[1]][0])
            Lo_counter.append(Lo_[idx[0]:idx[1]][0])

    MP['eps'] = c1*LT**2 * N2_overturn**1.5
    MP['N2_overturn'] = N2_overturn
    MP['LT'] = LT
    MP['Lo'] = Lo
    MP['LT_counter'] = LT_counter
    MP['Lo_counter'] = Lo_counter

    eps_tmean = np.nanmean(MP['eps'], axis=1)
    nans = np.isnan(eps_tmean)
    eps_tmean_ = np.zeros_like(eps_tmean)
    eps_tmean_[~nans] = eps_tmean[~nans]
    MP['eps_tzint'] = -1000*np.trapz(eps_tmean_, MP['z'][:, 0])

    nans = np.isnan(MP['eps'])
    eps_ = np.zeros_like(MP['eps'])
    eps_[~nans] = MP['eps'][~nans]
    MP['eps_zint'] = -1000*np.trapz(eps_, MP['z'][:, 0], axis=0)

    print('- Saving data')
    io.savemat(os.path.join(save_dir, name + '.mat'), MP)
