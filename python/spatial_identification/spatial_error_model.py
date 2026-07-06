"""
Spatial error identification for MCG spherical verification data.

Reconstructed reference implementation
---------------------------------------
The scripts that produced the original project's chapter-5 results
(resolution_numerique.py and the systeme_alpha / systeme_beta / systeme_gamma
helper functions referenced in the report) were not recovered with the
project archive - only screenshots of the results remain in the report.
This module is a clean re-implementation of the *documented* mathematical
model (report section 5.1-5.2), built directly from the design-matrix
formulation rather than from the original three-separate-single-parameter
routines.

The report explicitly flags that the original alpha/beta estimation was
inconsistent by one to two orders of magnitude (Table 5.1) and recommends
solving the full parameter vector jointly by least squares instead of
estimating each angle from an isolated one-parameter system. That is
exactly what this module does: a single joint fit of
    x = [u, v, w, dR, alpha, beta, gamma]
by SVD-based pseudoinverse, which is both simpler and avoids the
parameter-coupling issue that most likely caused the original discrepancy.

Model
-----
For a translated, radius-offset, slightly non-orthogonal sphere, the
scalar normal deviation returned by the CMM at point i is modeled to
first order as

    eps_i = u*nx_i + v*ny_i + w*nz_i + dR
            - alpha*z_i*ny_i - beta*x_i*nz_i - gamma*y_i*nx_i + eta_i

so that each row of the design matrix A is
    A_i = [nx_i, ny_i, nz_i, 1, -z_i*ny_i, -x_i*nz_i, -y_i*nx_i]

and the joint least-squares estimate is x_hat = pinv(A) @ eps.

Author: Dev Kumar, Thien Ho (model); reference re-implementation for
publication.
Project: Spatial Verification of a Coordinate Measuring Machine - ENSAM Lille
"""

import numpy as np

PARAM_NAMES = ['u', 'v', 'w', 'dR', 'alpha', 'beta', 'gamma']


def load_deviation_file(filename):
    """Load a CMM deviation file.

    Expected whitespace-separated columns, one point per row:
        x  y  z  nx  ny  nz  eps
    where (x, y, z) is the nominal point, (nx, ny, nz) the unit probing
    normal, and eps the scalar normal deviation reported by the CMM.

    Returns
    -------
    xyz : (N, 3) ndarray
    normals : (N, 3) ndarray
    eps : (N,) ndarray
    """
    data = np.loadtxt(filename)
    xyz = data[:, 0:3]
    normals = data[:, 3:6]
    eps = data[:, 6]
    return xyz, normals, eps


def build_design_matrix(xyz, normals, model='full'):
    """Build the least-squares design matrix A.

    Parameters
    ----------
    xyz : (N, 3) ndarray
        Nominal point coordinates.
    normals : (N, 3) ndarray
        Unit probing normals.
    model : {'translation_radius', 'full'}
        'translation_radius' fits x = [u, v, w, dR] only (section 5.1).
        'full' additionally fits the three small squareness terms
        alpha, beta, gamma (section 5.2).

    Returns
    -------
    A : (N, p) ndarray
    """
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    nx, ny, nz = normals[:, 0], normals[:, 1], normals[:, 2]
    ones = np.ones_like(nx)

    if model == 'translation_radius':
        return np.column_stack([nx, ny, nz, ones])
    elif model == 'full':
        return np.column_stack([
            nx, ny, nz, ones,
            -z * ny, -x * nz, -y * nx,
        ])
    else:
        raise ValueError("model must be 'translation_radius' or 'full'")


def fit(xyz, normals, eps, model='full'):
    """Solve the joint least-squares problem eps = A @ x by SVD.

    Returns a dict with the parameter estimate, residuals, residual
    variance, parameter covariance estimate, and the condition number of
    the design matrix (all quantities defined in report section 5.1).
    """
    A = build_design_matrix(xyz, normals, model=model)
    n, p = A.shape

    x_hat, residual_sum_sq, rank, singular_values = np.linalg.lstsq(A, eps, rcond=None)

    residuals = eps - A @ x_hat
    dof = max(n - p, 1)
    sigma2 = float(np.sum(residuals ** 2) / dof)

    # Cov(x_hat) ~= sigma^2 * (A^T A)^-1, computed via the pseudoinverse
    # of A so that no explicit matrix inversion is required.
    A_pinv = np.linalg.pinv(A)
    cov = sigma2 * (A_pinv @ A_pinv.T)

    condition_number = float(singular_values[0] / singular_values[-1]) if singular_values[-1] > 0 else np.inf

    names = PARAM_NAMES if model == 'full' else PARAM_NAMES[:4]

    return {
        'model': model,
        'names': names,
        'x_hat': dict(zip(names, x_hat)),
        'residuals': residuals,
        'residual_rms': float(np.sqrt(np.mean(residuals ** 2))),
        'sigma2': sigma2,
        'cov': cov,
        'std_errors': dict(zip(names, np.sqrt(np.diag(cov)))),
        'condition_number': condition_number,
        'singular_values': singular_values,
        'rank': rank,
    }


def generate_synthetic_dataset(radius, true_params, n_azimuth=60, n_rings=10,
                                elevation_deg=(-45, 45), noise_std=0.0, seed=0):
    """Generate a synthetic spherical dataset for a prescribed parameter
    vector, to validate `fit` end-to-end (compare estimated vs. prescribed
    parameters, as in report Table 5.1).

    Parameters
    ----------
    radius : float
        Nominal sphere radius, mm.
    true_params : dict
        Any subset of {'u','v','w','dR','alpha','beta','gamma'} in mm
        (u, v, w, dR) and radians (alpha, beta, gamma). Missing keys
        default to zero.
    noise_std : float
        Standard deviation of additive Gaussian noise on eps, mm.

    Returns
    -------
    xyz, normals, eps : as required by `fit`
    """
    rng = np.random.default_rng(seed)
    p = {k: true_params.get(k, 0.0) for k in PARAM_NAMES}

    phi = np.radians(np.linspace(elevation_deg[0], elevation_deg[1], n_rings))
    theta = 2 * np.pi * np.arange(n_azimuth) / n_azimuth

    PHI, THETA = np.meshgrid(phi, theta, indexing='ij')
    PHI = PHI.ravel()
    THETA = THETA.ravel()

    x = radius * np.cos(PHI) * np.cos(THETA)
    y = radius * np.cos(PHI) * np.sin(THETA)
    z = radius * np.sin(PHI)
    xyz = np.column_stack([x, y, z])

    normals = xyz / radius

    A_full = build_design_matrix(xyz, normals, model='full')
    x_true = np.array([p['u'], p['v'], p['w'], p['dR'], p['alpha'], p['beta'], p['gamma']])
    eps = A_full @ x_true
    if noise_std > 0:
        eps = eps + rng.normal(0, noise_std, size=eps.shape)

    return xyz, normals, eps


def _demo():
    """Reproduce a Table-5.1-style synthetic validation check."""
    true_params = {
        'u': 0.010, 'v': 0.020, 'w': 0.030, 'dR': 0.100,
        'alpha': 0.1e-6, 'beta': 5.0e-6, 'gamma': 10.0e-6,   # rad
    }

    xyz, normals, eps = generate_synthetic_dataset(radius=151.0, true_params=true_params)
    result = fit(xyz, normals, eps, model='full')

    print(f"{'Parameter':<10}{'Prescribed':>14}{'Estimated':>14}")
    for name in PARAM_NAMES:
        prescribed = true_params[name]
        estimated = result['x_hat'][name]
        if name in ('alpha', 'beta', 'gamma'):
            print(f"{name:<10}{prescribed*1e6:14.4f}{estimated*1e6:14.4f}   (microrad)")
        else:
            print(f"{name:<10}{prescribed:14.6f}{estimated:14.6f}   (mm)")

    print(f"\nCondition number kappa(A) = {result['condition_number']:.3e}")
    print(f"Residual RMS = {result['residual_rms']:.3e} mm (should be ~0 for noise-free data)")


if __name__ == '__main__':
    _demo()
