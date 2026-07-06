"""
Spherical-segment trajectory generation for Machine Checking Gauge (MCG)
verification of a coordinate measuring machine.

Generates M latitude rings evenly spaced in elevation between -45 deg
and +45 deg for each available MCG arm length, together with the
outward unit radial normal at every point, and exports a six-column
ASCII file (x, y, z, nx, ny, nz) per arm length for import into the CMM
measurement software.

Note on this implementation
----------------------------
The project's original point-generation script contained a commented-out
spherical-generation block parameterized by colatitude from the pole
(angle = j*(pi/20) + pi/4, j = 1..10). Working through that
parameterization shows it spans elevation angles from -45 deg up to only
+36 deg (an off-by-one asymmetry: 10 steps starting at j=1 rather than
j=0), not the symmetric +/-45 deg band described in the report.

This module instead implements the report's own, explicit, symmetric
elevation formula directly (section 3.3):
    phi_j = -45 deg + j * 90 deg / (M - 1),   j = 0 .. M-1
    z_j   = R * sin(phi_j)
    rho_j = R * cos(phi_j)
which is unambiguous and matches the stated +/-45 deg coverage.

Author: Dev Kumar, Thien Ho (original project and model);
        Python re-implementation for publication.
Project: Spatial Verification of a Coordinate Measuring Machine - ENSAM Lille
"""

import numpy as np
import matplotlib.pyplot as plt

ARM_LENGTHS_MM = [101, 151, 226, 380, 532, 685]

N_AZIMUTH = 60          # points per latitude ring
N_RINGS = 10            # latitude rings, M
ELEVATION_RANGE_DEG = (-45.0, 45.0)

CENTER = [0, 0, 0]


def generate_sphere(radius, center=CENTER, n_azimuth=N_AZIMUTH, n_rings=N_RINGS,
                     elevation_range_deg=ELEVATION_RANGE_DEG, uniform_density=False,
                     target_spacing=None):
    """Return the nominal points and outward unit normals of the
    spherical segment between elevation_range_deg[0] and [1].

    Parameters
    ----------
    radius : float
        Sphere radius (MCG arm length), mm.
    n_azimuth : int
        Azimuthal points per ring (used when uniform_density=False).
    n_rings : int
        Number of latitude rings, M, evenly spaced in elevation
        (both endpoints included).
    uniform_density : bool
        If True, scale the number of azimuthal points on each ring with
        its circumference (see report section 6.2) so that adjacent
        chord lengths remain approximately constant, instead of using a
        fixed n_azimuth on every ring.
    target_spacing : float, optional
        Target chord length (mm) when uniform_density=True. Defaults to
        the chord length implied by n_azimuth on the largest ring.

    Returns
    -------
    points : (N, 3) ndarray
    normals : (N, 3) ndarray
    """
    phi = np.radians(np.linspace(elevation_range_deg[0], elevation_range_deg[1], n_rings))

    if target_spacing is None:
        target_spacing = 2 * np.pi * radius / n_azimuth

    points = []
    normals = []

    for phi_j in phi:
        rho_j = radius * np.cos(phi_j)
        z_j = radius * np.sin(phi_j)

        n_j = n_azimuth
        if uniform_density:
            n_j = max(3, int(np.ceil(2 * np.pi * rho_j / target_spacing)))

        theta = 2 * np.pi * np.arange(1, n_j + 1) / n_j

        x = rho_j * np.cos(theta) + center[0]
        y = rho_j * np.sin(theta) + center[1]
        z = np.full(n_j, z_j + center[2])

        # Outward unit radial normal: (point - center) / radius
        nx = (x - center[0]) / radius
        ny = (y - center[1]) / radius
        nz = (z - center[2]) / radius

        points.append(np.column_stack([x, y, z]))
        normals.append(np.column_stack([nx, ny, nz]))

    return np.vstack(points), np.vstack(normals)


def export_points(filename, points, normals):
    """Write a six-column (x, y, z, nx, ny, nz) ASCII file for CALYPSO."""
    with open(filename, 'w') as f:
        for (x, y, z), (nx, ny, nz) in zip(points, normals):
            f.write(f'{x:.7f} {y:.7f} {z:.7f} {nx:.7f} {ny:.7f} {nz:.7f}\n')


def main():
    for radius in ARM_LENGTHS_MM:
        points, normals = generate_sphere(radius)
        export_points(f'coordinates_arm_{radius}mm_sphere.txt', points, normals)

    # 3-D preview of the largest-radius spherical segment
    radius = ARM_LENGTHS_MM[-1]
    points, normals = generate_sphere(radius)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(points[:, 0], points[:, 1], points[:, 2])
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_zlabel('Z (mm)')
    plt.show()


if __name__ == '__main__':
    main()
