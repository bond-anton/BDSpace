import numpy as np


def generate_cuboid(a=1, b=1, c=1, origin=np.array([0, 0, 0])):
    cube = np.array([[0, 0, 0], [0, b, 0], [a, 0, 0], [a, b, 0],
                     [0, 0, c], [0, b, c], [a, 0, c], [a, b, c]])
    return cube[:] - origin


def generate_sphere(resolution=180, radius=1):
    latitudes = np.linspace(0, np.pi, resolution * 0.5, endpoint=True)
    t = np.linspace(0, 2*np.pi, resolution, endpoint=True)
    points = np.ones((resolution * resolution * 0.5, 3), dtype=np.float)
    start = 0
    for latitude in latitudes:
        end = start + len(t)
        points[start:end, 2] *= radius * np.cos(latitude)
        points[start:end, 1] *= radius * np.sin(t) * np.sin(latitude)
        points[start:end, 0] *= radius * np.cos(t) * np.sin(latitude)
        start = end
    return points


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
