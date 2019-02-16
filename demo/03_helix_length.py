import numpy as np
from BDSpace.Coordinates import Cartesian
from BDSpace.Curve.Parametric import Helix
from matplotlib import pyplot as plt


def plot_tree(mesh_tree, ax=None):
    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'r', 'g', 'b', 'c', 'm', 'y', 'k',
              'r', 'g', 'b', 'c', 'm', 'y', 'k', 'r', 'g', 'b', 'c', 'm', 'y', 'k']
    styles = ['-', ':', '--', '-', ':', '--', '-', '-', ':', '--', '-', ':', '--', '-',
              '-', ':', '--', '-', ':', '--', '-', '-', ':', '--', '-', ':', '--', '-']
    show = False
    if ax is None:
        _, ax = plt.subplots(3, 1)
        show = True
    for level in mesh_tree.tree.keys():
        for i, mesh in enumerate(mesh_tree.tree[level]):
            ax[0].plot(mesh.physical_nodes, np.ones(mesh.num) * level, colors[level] + styles[i] + 'o')
            ax[1].plot(mesh.physical_nodes, mesh.solution, colors[level] + styles[i] + 'o')
            ax[2].plot(mesh.physical_nodes, mesh.residual, colors[level] + styles[i] + 'o')
    ax[0].set_ylim([-1, max(mesh_tree.tree.keys()) + 1])
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    if show:
        plt.show()


coordinate_system = Cartesian()
coordinate_system.rotate_axis_angle(np.ones(3, dtype=np.double), np.deg2rad(45))

left_helix = Helix(name='Left Helix', coordinate_system=coordinate_system,
                   radius=2, pitch=0.5, start=0, stop=np.pi * 2, right=False)

print('Helix length:', left_helix.length())
tree = left_helix.mesh_tree()

plot_tree(tree)
plt.show()

flat = tree.flatten()
print(np.sum(flat.solution))

step = 1e-10
t = np.linspace(0.0, 2 * np.pi, num=100)


xyz = np.asarray(left_helix.tangent(t))
plt.plot(t, xyz[:, 0], 'r-o')
plt.plot(t, xyz[:, 1], 'b-o')
plt.plot(t, xyz[:, 2], 'g-o')
plt.show()

dl = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
plt.plot(t, dl, 'b-o')
print(np.trapz(dl, t))
plt.show()
