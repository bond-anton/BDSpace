import numpy as np


def generate_cuboid(a=1, b=1, c=1, origin=np.array([0, 0, 0])):
    cube = np.array([[0, 0, 0], [0, b, 0], [a, 0, 0], [a, b, 0],
                     [0, 0, c], [0, b, c], [a, 0, c], [a, b, c]])
    return cube[:] - origin, (2, 2, 2)


def generate_sphere(r, theta, phi):
    """
    Generate points for structured grid for a spherical shell volume.
    This method is useful for generating a unstructured cylindrical mesh for VTK.
    :param r: radius grid
    :param phi: azimuthal angle grid
    :param theta: polar angle grid
    """
    # Find the x values and y values for each plane.
    x_plane = (np.sin(theta) * np.cos(phi) * r[:, None]).ravel()
    y_plane = (np.sin(theta) * np.cos(phi) * r[:, None]).ravel()
    z = (np.cos(theta) * r[:, None])

    # Allocate an array for all the points.  We'll have len(x_plane)
    # points on each plane, and we have a plane for each z value, so
    # we need len(x_plane)*len(z) points.
    points = np.empty([len(x_plane) * len(z), 3])

    # Loop through the points for each plane and fill them with the
    # correct x,y,z values.
    start = 0
    for z_plane in z:
        end = start+len(x_plane)
        # slice out a plane of the output points and fill it
        # with the x,y, and z values for this plane.  The x,y
        # values are the same for every plane.  The z value
        # is set to the current z
        plane_points = points[start:end]
        plane_points[:, 0] = x_plane
        plane_points[:, 1] = y_plane
        plane_points[:, 2] = z_plane
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
