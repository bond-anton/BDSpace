from __future__ import division
import numpy as np

from Space.Coordinates import Cartesian, transforms
from Space.Coordinates.transforms import unit_vector
from Space.Curve.Parametric import Helix, Arc


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
    arc_coordinate_system = Cartesian(basis=np.copy(coordinate_system.basis), origin=np.copy(global_point1),
                                      name='Arc coordinate_system')
    r_theta_phi = transforms.cartesian_to_spherical(direction)
    arc_coordinate_system.rotate_axis_angle([0, 0, 1], r_theta_phi[2])
    arc_coordinate_system.rotate_axis_angle([0, 1, 0], r_theta_phi[1])
    local_point2 = arc_coordinate_system.to_local(global_point2)
    print local_point2
    origin = np.array([np.sqrt((2 * radius)**2 - local_point2[2]**2), 0, 0])
    print origin
    direction = local_point2 - origin
    print direction
    print 'Arc coordinate_system origin offset:', origin
    distance = np.sqrt(np.dot(direction, direction))
    print distance
    origin = arc_coordinate_system.to_parent(origin)
    arc_coordinate_system.origin = np.copy(origin)
    r_theta_phi = transforms.cartesian_to_spherical(direction)
    arc_coordinate_system.rotate_axis_angle([0, 0, 1], r_theta_phi[2])
    arc_coordinate_system.rotate_axis_angle([0, 1, 0], r_theta_phi[2])
    local_point1 = arc_coordinate_system.to_local(coordinate_system.to_parent(point1)) - np.array([0, 0, radius])
    local_point2 = arc_coordinate_system.to_local(coordinate_system.to_parent(point2)) - np.array([0, 0, radius])
    print local_point2
    start = np.pi - np.arccos(local_point1[2] / radius)
    stop = np.pi - np.arccos(local_point2[2] / radius)
    print start, stop
    path = Arc(coordinate_system=arc_coordinate_system, radius=radius, start=start, stop=stop, right=right)
    return path
