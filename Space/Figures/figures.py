import numpy as np


def generate_cube_points(a=1, b=1, c=1, origin=np.array([0, 0, 0])):
    cube = np.array([[0, 0, 0], [0, b, 0], [a, 0, 0], [a, b, 0],
                     [0, 0, c], [0, b, c], [a, 0, c], [a, b, c]])
    return cube[:] - origin


def generate_cone(phi, z=np.array([0, 100.0]), theta=np.pi/4, hole=5, r_min=20):
    z_min = r_min / np.tan(theta)
    points = np.empty([2 * len(phi) * len(z), 3])
    start = 0
    for z_plane in z:
        r = np.array([hole, (z_min + z_plane) * np.tan(theta)])
        x_plane = (np.cos(phi) * r[:, None]).ravel()
        y_plane = (np.sin(phi) * r[:, None]).ravel()
        end = start + len(x_plane)
        plane_points = points[start:end]
        plane_points[:, 0] = x_plane
        plane_points[:, 1] = y_plane
        plane_points[:, 2] = z_plane
        start = end
    return points


def generate_cylinder(phi, z, r_outer, hole):
    points = np.empty([2 * len(phi) * len(z), 3])
    r = np.array([hole, r_outer])
    start = 0
    for z_plane in z:
        x_plane = (np.cos(phi) * r[:, None]).ravel()
        y_plane = (np.sin(phi) * r[:, None]).ravel()
        end = start + len(x_plane)
        plane_points = points[start:end]
        plane_points[:, 0] = x_plane
        plane_points[:, 1] = y_plane
        plane_points[:, 2] = z_plane
        start = end
    return points
