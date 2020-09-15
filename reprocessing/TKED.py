# -*- coding: utf-8 -*-


import numpy as np
import utils


def intermediate_profile1(x, hinge=1000, delta=1e-3, kind='bottom up'):
    """Generate an intermediate profile of some quantity. Ferron et. al. 1998.

    Parameters
    ----------
    x : 1D array
        Temperature or density.
    hinge : float, optional
        Hinge temperature or density.
    delta : float, optional
        Step,
    kind : string, optional
        Either 'bottom up', 'top down' or 'average'.

    Returns
    -------
    y : 1D array
        Reference buoyancy frequency [s-2]

    """
    xf = np.flipud(x)

    xtd = np.zeros_like(x)
    xbu = np.zeros_like(x)

    ntd = np.fix(x[0]/delta - hinge/delta)
    nbu = np.fix(xf[0]/delta - hinge/delta)

    xtd[0] = hinge + ntd*delta
    xbu[0] = hinge + nbu*delta

    for i in range(len(x) - 1):
        ntd = np.fix(x[i+1]/delta - xtd[i]/delta)
        nbu = np.fix(xf[i+1]/delta - xbu[i]/delta)

        xtd[i+1] = xtd[i] + ntd*delta
        xbu[i+1] = xbu[i] + nbu*delta

    xbu = np.flipud(xbu)

    xav = (xtd + xbu)/2.

    if 'up' in kind:
        return xbu
    elif 'down' in kind:
        return xtd
    elif 'av' in kind:
        return xav


def thorpe_scales1(z, x, acc=1e-3, R0=0.25, Nsq=None, full_output=False,
                   Nsq_method='bulk', use_int_prof=False, **ip_kwargs):
    """Estimate thorpe scales. Thorpe et. al. 1977
    Contains Gargett and Garner 2008 validation ratio.

    Parameters
    ----------
    z : 1D array
        Height. [m] (Negative depth!)
    x : 1D array
        Density. [kg m-3]
    acc : float, optional
        Accuracy of the x measurement.
    R0 : float, optional
        Validation ratio criteria, default 0.25.
    Nsq : 1D array, optional
        Buoyancy frequency squared. [rad2 s-2]
    full_output : boolean, optional
        Return all diagnostic variables. Also calculates N squared.
    Nsq_method : string, optional
        The method used to estimated buoyancy frequency. The options are
        'endpt' or 'bulk'. See Mater et. al. 2015 for a discussion.
    use_int_prof : boolean, optional
        Use the intermediate profile method of Ferron.
    ip_kwargs : dict, optional
        Keyword arguments for the intermediate profile method.

    Returns
    -------
    LT : 1D array
        Thorpe scales. [m]
    Td : 1D array
        Thorpe displacements. [m]
    Nsqo : 1D array, optional
        Buoyancy frequency of overturns. [rad2 s-2]
    Lo : 1D array, optional
        Overturn length. [m]
    R : 1D array, optional
        Overturn ratio.
    x_sorted : 1D array, optional
        Sorted density. [kg m-3]
    idxs : 1D array, optional
        Indexes required to sort.
    """
    g = -9.807  # Gravitational acceleration [m s-2]
    LT = np.zeros_like(x)
    Lo = np.zeros_like(x)
    R = np.zeros_like(x)
    Nsqo = np.zeros_like(x)

    # x should be increasing for this algorithm to work.
    flip_x = False
    if x[0] > x[-1]:
        x = np.flipud(x)
        z = np.flipud(z)
        if Nsq is not None:
            Nsq = np.flipud(Nsq)
        flip_x = True

    if use_int_prof:
        x = intermediate_profile(x, **ip_kwargs)

    # This is for estimating the length of the overturns.
    dz = 0.5*(z[:-2] - z[2:])
    dz = np.hstack((dz[0], dz, dz[-1]))

    # Make sure that no overturns involve the first or last points.
    x[0] = x.min() - 1e-4
    x[-1] = x.max() + 1e-4

    # Sort the profile.
    idxs = np.argsort(x)
    x_sorted = x[idxs]
    # Calculate thorpe displacements.
    Td = z[idxs] - z
    # Index displacements.
    idxs_disp = idxs - np.arange(len(idxs))
    # Overturn bounds where cumulative sum is zero.
    idxs_sum = np.cumsum(idxs_disp)
    # Find overturns.
    odxs_ = utils.contiguous_regions(idxs_sum > 0)

    if odxs_.size == 0:  # No oveturns at all!
        if full_output:
            return LT, Td, Nsqo, Lo, R, x_sorted, idxs
        else:
            return LT

    cut = (odxs_[:, 1] - odxs_[:, 0]) == 1
    if odxs_[0, 0] == 0:
        cut[0] = True

    odxs = odxs_[~cut, :]

    # Calculate the RMS thorpe displacement over each overturn.
    for j1, j2 in odxs:
        odx = slice(j1, j2)
        # Check for noise.
        q = x_sorted[j2] - x_sorted[j1]
        if q < acc:
            continue

        # Overturn ratio of Gargett & Garner
        Tdo = Td[odx]
        dzo = dz[odx]
        L_tot = np.sum(dzo)
        L_neg = np.sum(dzo[Tdo < 0])
        L_pos = np.sum(dzo[Tdo > 0])
        R_ = np.minimum(L_neg/L_tot, L_pos/L_tot)
        if R_ < R0:
            continue

        # Store data.
        Lo[odx] = L_tot
        R[odx] = R_
        LT_ = np.sqrt(np.mean(Tdo**2))
        LT[odx] = LT_
        if Nsq_method == 'endpt':
            Nsqo[odx] = -g*q/(np.mean(x_sorted[odx])*L_tot)
        elif Nsq_method == 'bulk':
            dx = x[odx] - x_sorted[odx]
            dxrms = np.sqrt(np.mean(dx**2))
            Nsqo[odx] = -g*dxrms/(np.mean(x_sorted[odx])*LT_)
        else:
            raise ValueError("Nsq_method can be either 'endpt' or 'bulk'.")

    # Lastly if the arrays were not increasing at the beginning and were
    # flipped they need to be put back how they were.
    if flip_x:
        LT = np.flipud(LT)
        Td = np.flipud(Td)
        x_sorted = np.flipud(x_sorted)
        idxs = np.flipud(idxs)
        Lo = np.flipud(Lo)
        R = np.flipud(R)
        Nsqo = np.flipud(Nsqo)

    if full_output:
        return LT, Td, Nsqo, Lo, R, x_sorted, idxs
    else:
        return LT
