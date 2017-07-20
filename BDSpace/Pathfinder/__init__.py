from __future__ import division, print_function
import numpy as np

from BDSpace.Coordinates import Cartesian, transforms
from BDSpace.Coordinates.transforms import unit_vector
from BDSpace.Curve.Parametric import Line, Helix, Arc


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
    origin = coordinate_system.to_parent(point1)
    helix_coordinate_system = Cartesian(basis=np.copy(coordinate_system.basis), origin=np.copy(origin),
                                        name='Helix coordinate system')
    r_theta_phi = transforms.cartesian_to_spherical(direction)
    helix_coordinate_system.rotate_axis_angle([0, 0, 1], r_theta_phi[2])
    helix_coordinate_system.rotate_axis_angle([0, 1, 0], r_theta_phi[1])
    pitch = distance / int(loops)
    name = 'Right Helix' if right else 'Left Helix'
    path = Helix(name=name, coordinate_system=helix_coordinate_system,
                 radius=radius, pitch=pitch, start=0, stop=np.pi * 2 * int(loops), right=right)
    return path


def arc_between_two_points(coordinate_system, point1, point2, radius=1, right=True):
    global_point1 = coordinate_system.to_parent(point1)
    global_point2 = coordinate_system.to_parent(point2)
    direction = point2 - point1
    distance = np.sqrt(np.dot(direction, direction))
    arc_coordinate_system = Cartesian(basis=np.copy(coordinate_system.basis), origin=np.copy(global_point1),
                                      name='Arc coordinate_system')

    r_theta_phi = transforms.cartesian_to_spherical(direction)
    arc_coordinate_system.rotate_axis_angle([0, 0, 1], r_theta_phi[2])
    arc_coordinate_system.rotate_axis_angle([0, 1, 0], r_theta_phi[1] + np.pi/2)
    x_offset = -distance / 2
    y_offset = np.sqrt(radius**2 - x_offset**2)
    if right:
        y_offset *= -1
    arc_coordinate_system.origin = arc_coordinate_system.to_parent([x_offset, y_offset, 0])
    local_point1 = arc_coordinate_system.to_local(global_point1)
    local_point2 = arc_coordinate_system.to_local(global_point2)
    start = transforms.cartesian_to_spherical(local_point1)[2]
    stop = transforms.cartesian_to_spherical(local_point2)[2]
    if not right:
        start = 2 * np.pi - start
        stop = 2 * np.pi - stop
    path = Arc(coordinate_system=arc_coordinate_system, a=radius, b=radius, start=start, stop=stop, right=right)
    return path
