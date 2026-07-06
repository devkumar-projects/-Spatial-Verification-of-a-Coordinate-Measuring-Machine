"""
Diagnostic smoothing filters for exploratory analysis of CMM deviation
data (report chapter 6.4).

Reconstructed reference implementation
---------------------------------------
The original exploratory MATLAB script (Lissage.m) was not recovered
with the project archive; only a screenshot remains in the report
(figure 6.2), which applies a moving-average filter to the coordinate
differences between a theoretical and an experimental point file. This
module reproduces that same behaviour in Python, plus the first-order
exponential (IIR) filter given analytically in report section 6.4,
    y_f[k] = lambda * y_f[k-1] + (1 - lambda) * y[k],   lambda = exp(-Ts/tau)

IMPORTANT - filtering and conformity (report section 6.4)
------------------------------------------------------------
Smoothing can reduce peaks and alter circularity or sphericity results.
Raw data must always be retained. Unfiltered results should be used for
conformity decisions unless the filter is explicitly defined by the
applicable specification and its influence is included in the
uncertainty budget. These filters are appropriate for diagnostic
visualization and frequency-content analysis only.

Author: Dev Kumar, Thien Ho (original MATLAB routine); Python
reference re-implementation for publication.
Project: Spatial Verification of a Coordinate Measuring Machine - ENSAM Lille
"""

import numpy as np


def moving_average(signal, window_size=5):
    """Centered moving-average filter, equivalent to MATLAB's movmean.

    Parameters
    ----------
    signal : (N,) array-like
    window_size : int

    Returns
    -------
    smoothed : (N,) ndarray
    """
    signal = np.asarray(signal, dtype=float)
    kernel = np.ones(window_size) / window_size
    # 'same' mode with edge padding reproduces movmean's shrinking window
    # behaviour near the boundaries closely enough for diagnostic use.
    padded = np.pad(signal, (window_size // 2, window_size // 2), mode='edge')
    smoothed = np.convolve(padded, kernel, mode='valid')
    return smoothed[:len(signal)]


def exponential_filter(signal, tau, sample_time):
    """First-order (IIR) exponential smoothing filter.

    Implements y_f[k] = lambda*y_f[k-1] + (1-lambda)*y[k] with
    lambda = exp(-Ts/tau), matching the continuous-time transfer
    function H(s) = 1 / (1 + tau*s) (report section 6.4).

    Parameters
    ----------
    signal : (N,) array-like
    tau : float
        Filter time constant (same time unit as sample_time).
    sample_time : float
        Sampling interval Ts.

    Returns
    -------
    filtered : (N,) ndarray
    """
    signal = np.asarray(signal, dtype=float)
    lam = np.exp(-sample_time / tau)

    filtered = np.empty_like(signal)
    filtered[0] = signal[0]
    for k in range(1, len(signal)):
        filtered[k] = lam * filtered[k - 1] + (1 - lam) * signal[k]
    return filtered


def compare_theoretical_experimental(theoretical_file, experimental_file,
                                      window_size=5, output_file=None):
    """Reproduce the original Lissage.m workflow: load a theoretical and
    an experimental six-column (x,y,z,nx,ny,nz) point file, compute the
    coordinate-wise differences, smooth each difference channel with a
    moving average, and optionally write the corrected experimental
    coordinates to output_file.

    Returns
    -------
    corrected : (N, 6) ndarray
        Experimental coordinates with the smoothed difference removed.
    """
    theoretical = np.loadtxt(theoretical_file)
    experimental = np.loadtxt(experimental_file)

    diff = experimental - theoretical
    diff_smooth = np.column_stack([
        moving_average(diff[:, col], window_size) for col in range(diff.shape[1])
    ])

    corrected = experimental - diff_smooth

    if output_file:
        np.savetxt(output_file, corrected, fmt='%.7f')

    return corrected
