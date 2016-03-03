from __future__ import division
import numpy as np

from Space.Coordinates import Cartesian
from Space.Coordinates.transforms import unit_vector
from Space.Curve.Parametric import Helix, Arc


def helix_between_two_points(coordinate_system, point1, point2, r=1, loops=1, right=True):
    distance = np.sqrt(np.dot(point2 - point1, point2 - point1))
    direction = unit_vector(point2 - point1)
    print direction
    origin = coordinate_system.to_parent(point1)
    helix_coordinate_system = Cartesian(basis=np.copy(coordinate_system.basis), origin=np.copy(origin),
                                        name='Helix coordinate system')
    R = np.sqrt(direction[0]**2 + direction[1]**2)
    theta = np.arcsin(R)
    if R == 0:
        phi = 0
    else:
        phi = np.arccos(direction[0] / R)
    print R, np.rad2deg(theta), np.rad2deg(phi)
    helix_coordinate_system.rotate_axis_angle(helix_coordinate_system.basis[:, 2], -phi)
    helix_coordinate_system.rotate_axis_angle(helix_coordinate_system.basis[:, 1], -theta)
    h = distance / int(loops)
    path = Helix(name='Right Helix', coordinate_system=helix_coordinate_system,
                 r=r, h=h, start=0, stop=np.pi * 2 * int(loops), right=right)
    return path


def arc_between_two_points(coordinate_system, point1, point2, r=1, right=True):
    direction = unit_vector(point2 - point1)
    origin = coordinate_system.to_parent(point1)
    arc_coordinate_system = Cartesian(basis=np.copy(coordinate_system.basis), origin=np.copy(origin),
                                      name='Arc coordinate_system')
    R = np.sqrt(direction[0]**2 + direction[1]**2)
    theta = np.arcsin(R)
    if R == 0:
        phi = 0
    else:
        phi = np.arccos(direction[0] / R)
    arc_coordinate_system.rotate_axis_angle(arc_coordinate_system.basis[:, 2], -phi)
    arc_coordinate_system.rotate_axis_angle(arc_coordinate_system.basis[:, 1], -theta)
    local_point2 = arc_coordinate_system.to_local(coordinate_system.to_parent(point2)).ravel()
    origin2 = np.array([np.sqrt((2 * r)**2 - local_point2[2]**2), 0, 0])
    print 'Arc coordinate_system origin offset:', origin2
    distance = np.sqrt(np.dot(local_point2 - origin2, local_point2 - origin2))
    print distance
    direction = unit_vector(local_point2 - origin2)
    print direction
    origin = arc_coordinate_system.to_parent(origin2)
    arc_coordinate_system = Cartesian(basis=np.copy(arc_coordinate_system.basis), origin=np.copy(origin),
                                      name='Arc coordinate_system')
    R = np.sqrt(direction[0]**2 + direction[1]**2)
    theta = np.arcsin(R)
    if R == 0:
        phi = 0
    else:
        phi = np.arccos(direction[0] / R)
    print R, np.rad2deg(theta), np.rad2deg(phi)
    arc_coordinate_system.rotate_axis_angle(arc_coordinate_system.basis[:, 2], -phi)
    arc_coordinate_system.rotate_axis_angle(arc_coordinate_system.basis[:, 1], -theta)
    local_point1 = arc_coordinate_system.to_local(coordinate_system.to_parent(point1)).ravel() - np.array([0, 0, r]).ravel()
    local_point2 = arc_coordinate_system.to_local(coordinate_system.to_parent(point2)).ravel() - np.array([0, 0, r]).ravel()
    start = np.pi - np.arccos(local_point1[2] / r)
    stop = np.pi - np.arccos(local_point2[2] / r)
    print start, stop
    path = Arc(coordinate_system=arc_coordinate_system, r=r, start=start, stop=stop, right=True)
    return path
