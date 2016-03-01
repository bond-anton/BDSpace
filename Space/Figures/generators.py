import numpy as np


def generate_cuboid(a=1, b=1, c=1, origin=np.array([0, 0, 0])):
    """
    generates cuboid with origin offset
    :param a: length parameter a
    :param b: length parameter b
    :param c: length parameter c
    :param origin: offset point of cuboid center
    :return: grid points, dimensions tuple
    """
    cube = np.array([[0, 0, 0], [0, b, 0], [a, 0, 0], [a, b, 0],
                     [0, 0, c], [0, b, c], [a, 0, c], [a, b, c]])
    dims = (2, 2, 2)
    return cube[:] - origin, dims


def generate_sphere(phi, theta, r):
    """
    Generate points for structured grid for a spherical shell volume.
    This method is useful for generating a structured cylindrical mesh for VTK.
    :param phi: azimuthal angle array
    :param theta: polar angle array
    :param r: radius of sphere
    :return: grid points, dimensions tuple
    """
    r = np.array(r, dtype=np.float)
    points = np.empty([len(phi) * len(r) * len(theta), 3])
    start = 0
    for th in theta:
        x_plane = (np.cos(phi) * r[:, None] * np.sin(th)).ravel()
        y_plane = (np.sin(phi) * r[:, None] * np.sin(th)).ravel()
        z_plane = (np.ones_like(phi) * r[:, None] * np.cos(th)).ravel()
        end = start + len(x_plane)
        plane_points = points[start:end]
        plane_points[:, 0] = x_plane
        plane_points[:, 1] = y_plane
        plane_points[:, 2] = z_plane
        start = end
    dims = (len(phi), len(r), len(theta))
    return points, dims


def generate_cylinder(phi, z, r):
    """
    Generate points for structured grid for a cylindrical shell volume.
    This method is useful for generating a structured cylindrical mesh for VTK.
    :param phi: azimuthal angle array
    :param z: cylinder height array
    :param r: radius of cylinder
    :return: grid points, dimensions tuple
    """
    r = np.array(r, dtype=np.float)
    points = np.empty([len(phi) * len(r) * len(z), 3])
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
    dims = (len(phi), len(r), len(z))
    return points, dims


def generate_cone(phi, z=np.array([0, 100.0]), theta=np.pi/4, z_offset=5, r_min=20):
    """
    Generates structured grid of truncated hollow thick cone
    :param phi: azimuthal angles array
    :param z: height array
    :param theta: apex angle
    :param z_offset: thickness of cone wall
    :param r_min: truncation radius
    :return: grid points, dimensions tuple
    """
    z_min = r_min / np.tan(theta)
    if z_offset == 0:
        points = np.empty([1 * len(phi) * len(z), 3])
    else:
        points = np.empty([2 * len(phi) * len(z), 3])
    start = 0
    for z_plane in z:
        if z_offset != 0:
            r = np.array([(z_min + z_plane + z_offset) * np.tan(theta), (z_min + z_plane) * np.tan(theta)])
        else:
            r = np.array([(z_min + z_plane) * np.tan(theta)])
        x_plane = (np.cos(phi) * r[:, None]).ravel()
        y_plane = (np.sin(phi) * r[:, None]).ravel()
        end = start + len(x_plane)
        plane_points = points[start:end]
        plane_points[:, 0] = x_plane
        plane_points[:, 1] = y_plane
        plane_points[:, 2] = z_plane
        start = end
    dims = (len(phi), len(r), len(z))
    return points, dims


def generate_cone_with_cylindrical_hole(phi, z=np.array([0, 100.0]), theta=np.pi/4, hole_radius=5, r_min=20):
    """
    Generates structured grid of truncated cone with a throughout cylindrical hole along main axis
    :param phi: azimuthal angles array
    :param z: height array
    :param theta: apex angle
    :param hole_radius: radius of a hole
    :param r_min: truncation radius
    :return: grid points, dimensions tuple
    """
    if r_min < hole_radius:
        r_min = hole_radius
    z_min = r_min / np.tan(theta)
    points = np.empty([2 * len(phi) * len(z), 3])
    start = 0
    for z_plane in z:
        r = np.array([hole_radius, (z_min + z_plane) * np.tan(theta)])
        x_plane = (np.cos(phi) * r[:, None]).ravel()
        y_plane = (np.sin(phi) * r[:, None]).ravel()
        end = start + len(x_plane)
        plane_points = points[start:end]
        plane_points[:, 0] = x_plane
        plane_points[:, 1] = y_plane
        plane_points[:, 2] = z_plane
        start = end
    dims = (len(phi), len(r), len(z))
    return points, dims
