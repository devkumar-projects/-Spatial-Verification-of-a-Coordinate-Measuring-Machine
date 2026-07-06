"""
Planar (circular) trajectory generation for Machine Checking Gauge (MCG)
verification of a coordinate measuring machine.

For each available MCG arm length, generates N evenly spaced points on a
circle in the XY plane together with the outward unit probing normal at
each point, and exports them as a six-column ASCII file
(x, y, z, nx, ny, nz) that can be imported directly into the CMM
measurement software (ZEISS CALYPSO) as a nominal-point / normal-vector
program.

Point-count rationale
----------------------
Between two consecutive angular stations separated by dtheta = 2*pi/N,
the CMM moves along the chord rather than along the circular arc. For
the stylus to stay within the MCG's retaining fork, the chord length
    d_c = 2*R*sin(dtheta/2)
must remain below the fork's maximum admissible travel. A scan over
N = 50..59 at the longest arm (R = 685 mm, worst case) showed that
N = 30 is already sufficient; N = 60 was retained as a safety margin.

Author: Dev Kumar, Thien Ho
Project: Spatial Verification of a Coordinate Measuring Machine - ENSAM Lille
"""

import numpy as np
import matplotlib.pyplot as plt

# Available MCG nominal arm lengths, mm
ARM_LENGTHS_MM = [101, 151, 226, 380, 532, 685]

# Available MCG support heights, mm (not used in the planar campaign)
SUPPORT_HEIGHTS_MM = [31.75, 72, 127, 235]

N_POINTS = 60                      # points per revolution (see rationale above)
THETA_STEP = 2 * np.pi / N_POINTS  # angle between two consecutive points, rad

CENTER = [0, 0, 0]   # nominal X, Y, Z of the circle center
SUPPORT_HEIGHT = 0   # selected support height, mm


def generate_circle(radius, center=CENTER, n_points=N_POINTS, clockwise=True):
    """Return the nominal points and outward unit normals of a circle.

    Parameters
    ----------
    radius : float
        Circle radius (MCG arm length), mm.
    center : sequence of 3 floats
        Circle center coordinates (X, Y, Z).
    n_points : int
        Number of points on the circle.
    clockwise : bool
        Direction of traversal. The sign of the Y component is flipped
        to produce the complementary (counter-clockwise) path, which is
        used together with the clockwise path for the direction-history
        (hysteresis) comparison described in the report.

    Returns
    -------
    points : (n_points, 3) ndarray
    normals : (n_points, 3) ndarray
    """
    sign = 1 if clockwise else -1
    theta = THETA_STEP * np.arange(1, n_points + 1)

    x = radius * np.cos(theta) + center[0]
    y = sign * (-radius * np.sin(theta)) + center[1]
    z = np.full(n_points, center[2] + SUPPORT_HEIGHT)

    nx = np.cos(theta) + center[0]
    ny = sign * (-np.sin(theta)) + center[1]
    nz = np.full(n_points, center[2])

    points = np.column_stack([x, y, z])
    normals = np.column_stack([nx, ny, nz])
    return points, normals


def export_points(filename, points, normals):
    """Write a six-column (x, y, z, nx, ny, nz) ASCII file for CALYPSO."""
    with open(filename, 'w') as f:
        for (x, y, z), (nx, ny, nz) in zip(points, normals):
            f.write(f'{x:.7f} {y:.7f} {z:.7f} {nx:.7f} {ny:.7f} {nz:.7f}\n')


def main():
    for radius in ARM_LENGTHS_MM:
        points, normals = generate_circle(radius)
        export_points(f'coordinates_arm_{radius}mm_plane_normal.txt', points, normals)

    # Preview the largest-radius circle with its outward normals
    radius = ARM_LENGTHS_MM[-1]
    points, normals = generate_circle(radius)

    plt.scatter(points[:, 0], points[:, 1], label='Nominal circle')
    plt.quiver(points[:, 0], points[:, 1], normals[:, 0], normals[:, 1], color='r')
    plt.axis('equal')
    plt.xlabel('X (mm)')
    plt.ylabel('Y (mm)')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
