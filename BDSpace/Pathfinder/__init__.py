import numpy as np

from BDSpace.Coordinates import Cartesian, transforms
from BDSpace.Coordinates.transforms import unit_vector
from BDSpace.Curve import Line, Arc, Helix


def line_between_two_points(coordinate_system, point1, point2):
    direction = point2 - point1
    distance = np.sqrt(np.dot(direction, direction))
    v = unit_vector(direction)
    line_coordinate_system = Cartesian(basis=np.copy(coordinate_system.basis), origin=np.copy(coordinate_system.origin),
                                       name='Line path coordinate system')
    path = Line(name='Line Path', coordinate_system=line_coordinate_system,
                origin=point1, a=v[0], b=v[1], c=v[2],
                start=0, stop=distance)
    return path


def helix_between_two_points(coordinate_system, point1, point2, radius=1, loops=1, right=True):
    direction = point2 - point1
    distance = np.sqrt(np.dot(direction, direction))
    origin = coordinate_system.to_parent_vector(point1)
    helix_coordinate_system = Cartesian(basis=np.copy(coordinate_system.basis), origin=np.copy(origin),
                                        name='Helix coordinate system')
    r_theta_phi = transforms.cartesian_to_spherical_point(direction)
    helix_coordinate_system.rotate_axis_angle(np.array([0, 0, 1], dtype=np.double), r_theta_phi[2])
    helix_coordinate_system.rotate_axis_angle(np.array([0, 1, 0], dtype=np.double), r_theta_phi[1])
    pitch = distance / int(loops)
    name = 'Right Helix' if right else 'Left Helix'
    path = Helix(name=name, coordinate_system=helix_coordinate_system,
                 radius=radius, pitch=pitch, start=0, stop=np.pi * 2 * int(loops), right=right)
    return path


def arc_between_two_points(coordinate_system, point1, point2, radius=1, right=True):
    global_point = coordinate_system.to_parent(np.vstack((point1, point2)))
    direction = point2 - point1
    distance = np.sqrt(np.dot(direction, direction))
    arc_coordinate_system = Cartesian(basis=np.copy(coordinate_system.basis), origin=np.copy(global_point[0]),
                                      name='Arc coordinate_system')

    r_theta_phi = transforms.cartesian_to_spherical_point(direction)
    arc_coordinate_system.rotate_axis_angle(np.array([0, 0, 1], dtype=np.double), r_theta_phi[2])
    arc_coordinate_system.rotate_axis_angle(np.array([0, 1, 0], dtype=np.double), r_theta_phi[1] + np.pi/2)
    x_offset = -distance / 2
    y_offset = np.sqrt(radius**2 - x_offset**2)
    if right:
        y_offset *= -1
    arc_coordinate_system.origin = arc_coordinate_system.to_parent_vector(np.array([x_offset, y_offset, 0],
                                                                                   dtype=np.double))
    local_point = arc_coordinate_system.to_local(global_point)
    phi = transforms.cartesian_to_spherical(local_point)[:, 2]
    if right:
        start = phi[0]
        stop = phi[1]
    else:
        start = 2 * np.pi - phi[0]
        stop = 2 * np.pi - phi[1]
    path = Arc(coordinate_system=arc_coordinate_system, a=radius, b=radius, start=start, stop=stop, right=right)
    return path
