# CODE FOR ADVANCE TOPICS, SOLVES A PARABOLIC PARTIAL DIFFERENTIAL EQUATION ON A 2D SQUARE DOMAIN

# IMPORT OF LIBRARIES
import matplotlib.pyplot as plt
import numpy as np  # arrays and numerical application
from mpi4py import MPI
from scipy import sparse  # sparse matrices
from scipy.sparse import linalg
import math
from matplotlib import cm
from matplotlib.ticker import LinearLocator

# RESOURCE USAGE
import time
import tracemalloc


# FUNCTIONS


def create_simple_mesh(x, y):
    """Function used to create a simple mesh of a squared domain using the coordinates x and y, it also
    creates the s vector for each element to be used for the material specific properties"""
    side = len(x)
    xv = np.tile(x, side)  # X coordinates in row by row
    yv = np.repeat(y, side)  # Y coordinates in column*row
    nodes = np.asarray([xv, yv]).transpose()

    topology = np.zeros(shape=(((side - 1) ** 2) * 2, 3))
    ind = 0
    nd = 1

    while nd <= len(nodes):  # cycle over nodes
        if nd % side != 0 and nd <= side * (side - 1):  # skip nodes on the right edge
            topology[ind, :] = [nd, nd + side + 1, nd + side]  # top of upper triang
            topology[ind + 1, :] = [nd, nd + 1, nd + side + 1]  # top of lower triang
            ind = ind + 2  # two rows at a time
        nd = nd + 1  # go to next node

    topology = np.array(topology.astype(int))

    return nodes, topology


def const_voronoi_cells(top, area, nds):
    """Calculates the voronoi cells associated to each element assuming triangular regular mesh with
    all elements being of area area"""
    vor_cells = np.zeros(len(nds))
    for i in range(len(nds)):
        vor_cells[i] = np.count_nonzero(top == i + 1) * area / 3

    return vor_cells


def hf_calc(nds, top, area, vor_cells, H):
    """Function used to assembly the P and H matrices element wise using the elements in
    top and the nodes in nds. The parameter sv is the vector containing the S values for each element"""

    fvt = np.zeros(len(nds))
    # This will be a point to split, how many processes can I use? How should I split?
    for nel, elem in enumerate(top):

        h_loc = local_H(nds, elem, area)
        f_loc = local_f(elem, vor_cells, nds)
        for i in range(3):
            ind_i = elem[i] - 1

            fvt[ind_i] = fvt[ind_i] + f_loc[i]

            for j in range(3):
                ind_j = elem[j] - 1
                H[ind_i, ind_j] = H[ind_i, ind_j] + h_loc[i, j]

    return H, fvt


def local_H(nds, elem, area):
    """Function used to calculate the local stiffness matrix of an element"""

    ni = elem[0]
    nj = elem[1]
    nk = elem[2]

    b = np.array([nds[nj - 1, 1] - nds[nk - 1, 1], nds[nk - 1, 1] - nds[ni - 1, 0], nds[ni - 1, 1] - nds[nj - 1, 1]])
    c = np.array([nds[nk - 1, 0] - nds[nj - 1, 0], nds[ni - 1, 0] - nds[nk - 1, 0], nds[nj - 1, 0] - nds[ni - 1, 0]])
    loc_h = np.zeros(shape=(3, 3))

    for i in range(3):
        for j in range(3):
            loc_h[i, j] = (b[i] * b[j] + c[i] * c[j]) / (4 * area)

    return loc_h


def local_f(elem, vor_cells, nds):
    """Function to calculate the actual value of f element-wise"""

    ni = elem[0]
    nj = elem[1]
    nk = elem[2]

    loc_f = np.zeros(shape=(3, 1))

    for i, ni in enumerate([ni, nj, nk]):
        loc_f[i] = f_fun(nds[ni - 1, :]) * vor_cells[ni - 1]
    return loc_f


def dirichlet_bc(H, f, uv, dbc_nodes, penalty=10 ** 20):
    """Function to add the Dirichlet BCs to the problem: H is the sparse stiffness matrix, f the vector of known values,
    xv and yv the x and y coordinates of the nodes respectively and u the values of the solution"""

    for node in dir_nodes:
        H[node, node] = penalty
        f[node] = uv[node] * penalty

    return f, H


def u_fun(node):
    """Function to calculate the true value of the unknown function on the node node needed for the boundary
    conditions """
    x = node[0]
    y = node[1]

    return x ** 2 + y ** 2 - (x ** 2) * (y ** 2) - 1  # np.cos(2 * x ** 2 - y) + y ** 2


def f_fun(node):
    """Function to calculate the true value of the unknown function on the node node needed for the boundary
    conditions """
    x = node[0]
    y = node[1]

    return 4 - 2 * x ** 2 - 2 * y ** 2  # -(-4 * np.sin(2 * x ** 2 - y) - (16 * x ** 2 + 1) * np.cos(2 * x ** 2 - y) + 2)


# MPI INIZIALIZATION
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
root_pr = 0
assert size > 1

# TIME AND RESOURCES
if rank == root_pr:
    tracemalloc.start()
    starttime = time.time()

# DEFINE THINGS NOT CALCULATED EVERYWHERE
sendnodes = None
sendtop = None
uv_true = None
nodes = None
area = None
vor_cells = None
fv = None

# DOMAIN INIZIALIZATION
extr = 10
step = 0.1

# SOME EARLY CALCULATIONS ON THE PROBLEM SIZE
# NODES
side = 2 * int(extr / step) + 1
numnodes = side ** 2
nodes_per_process = math.floor(numnodes / size)
extra_nodes = math.ceil(numnodes / size)
countnodes = np.repeat(2 * nodes_per_process, size - 1)
countnodes = np.append(np.array(2 * extra_nodes), countnodes)
dispnodes = np.arange(2 * extra_nodes, 2 * numnodes, 2 * nodes_per_process)
dispnodes = np.append(0, dispnodes)

# SOLUTION
countsol = np.repeat(nodes_per_process, size - 1)
countsol = np.append(np.array(extra_nodes), countsol)
dispsol = np.arange(extra_nodes, numnodes, nodes_per_process)
dispsol = np.append(0, dispsol)

# TOPOLOGY
numelem = ((side - 1) ** 2) * 2
elem_per_process = math.floor(numelem / size)
extra_elem = math.ceil(numelem / size)
countelem = np.repeat(3 * elem_per_process, size - 1)
countelem = np.append(np.array(3 * extra_elem), countelem)
dispelem = np.arange(3 * extra_elem, 3 * numelem, 3 * elem_per_process)
dispelem = np.append(0, dispelem)

# VORONOI CELLS
countvor = np.repeat(elem_per_process, size - 1)
countvor = np.append(np.array(extra_elem), countvor)
dispvor = np.arange(extra_elem, numelem, elem_per_process)
dispvor = np.append(0, dispvor)

if rank == root_pr:
    x = np.arange(-extr, extr + step, step)
    y = np.arange(-extr, extr + step, step)

    nodes, topology = create_simple_mesh(x, y)

    shape = np.array([[1, nodes[topology[0, 0] - 1, 0], nodes[topology[0, 0] - 1, 1]],
                      [1, nodes[topology[0, 1] - 1, 0], nodes[topology[0, 1] - 1, 1]],
                      [1, nodes[topology[0, 2] - 1, 0], nodes[topology[0, 2] - 1, 1]]])

    area = np.linalg.det(shape) / 2
    vor_cells = const_voronoi_cells(topology, area, nodes)
    sendnodes = nodes.flatten().astype("float64")
    recvnodes = np.empty([extra_nodes, 2], dtype=np.float64)
    uv_true = np.empty(numnodes, dtype=np.float64)
    fv = np.empty(numnodes, dtype=np.float64)
else:
    recvnodes = np.empty([nodes_per_process, 2], dtype=np.float64)

comm.Scatterv([sendnodes, countnodes, dispnodes, MPI.DOUBLE], recvnodes, root=root_pr)
locnodes = recvnodes
uv_true_loc = np.apply_along_axis(u_fun, 1, locnodes)

comm.Gatherv(uv_true_loc, [uv_true, countsol, dispsol, MPI.DOUBLE], root=root_pr)

# SENDING DATA NEEDED FOR CALCULATING H AND fv
nodes = comm.bcast(nodes, root=root_pr)
area = comm.bcast(area, root=root_pr)
vor_cells = comm.bcast(vor_cells, root=root_pr)
H = sparse.lil_matrix((numnodes, numnodes))

if rank == root_pr:
    sendtop = topology.flatten().astype('int32')
    recvtop = np.empty([extra_elem, 3], dtype=np.int32)
else:
    recvtop = np.empty([elem_per_process, 3], dtype=np.int32)

comm.Scatterv([sendtop, countelem, dispelem, MPI.INT], recvtop, root=root_pr)
loc_topology = recvtop

H = comm.bcast(H, root=root_pr)
H_loc, fv_loc = hf_calc(nodes, loc_topology, area, vor_cells, H)

if rank == root_pr:
    H = H_loc.tocsr()
    fv = fv_loc

for r in range(1, size):
    if rank == r:
        comm.send(H_loc, dest=0, tag=r)
for r in range(1, size):
    if rank == root_pr:
        H_loc = comm.recv(source=r, tag=r)
        H = H + H_loc.tocsr()

for r in range(1, size):
    if rank == r:
        comm.send(fv_loc, dest=0, tag=r)
for r in range(1, size):
    if rank == root_pr:
        fv_loc = comm.recv(source=r, tag=r)
        fv = fv + fv_loc


if rank == root_pr:
    dir_nodes = np.array(np.concatenate((range(side), range(side * (side - 1), side ** 2, 1),
                                         range(side, side * (side - 1), side),
                                         range(2 * side - 1, side * (side - 1), side))))
    # BC plot
    plt.scatter(nodes[:, 0], nodes[:, 1])
    plt.scatter(nodes[dir_nodes[:], 0], nodes[dir_nodes[:], 1])
    plt.title('NODES OF THE FE MESH')
    plt.xlabel('x axis')
    plt.ylabel('y axis')
    plt.legend(['Computed nodes', 'Dirichlet BC nodes'])
    plt.savefig('bcplot.png')

    fv, H = dirichlet_bc(H, fv, uv_true, dir_nodes, penalty=10 ** 15)

    uv, info = linalg.cg(H, fv, atol=10 ** (-50), tol=10 ** (-50))

    maxerr = np.max(np.abs(uv - uv_true))

    # CODE FOR PLOT
    x = np.array([x] * len(x))
    y = np.array([y] * len(x)).transpose()
    z = uv.reshape(len(x), len(x))

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    # Plot the surface.
    surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')

    plt.savefig('Solution of the PDE')

    endtime = time.time()
    currentmem, peakmem = tracemalloc.get_traced_memory()
    f = open("results.txt", "w")
    f.write('RESULTS:\n')
    f.write('max absolute error: ' + str(maxerr) + '\n')
    f.write('Current used memory: ' + str(currentmem / (1024 ** 2)) + ' MB\nPeak used memory: ' + str(
        peakmem / (1024 ** 2)) + ' MB\n')
    f.write('Total elapsed time: ' + str(endtime - starttime) + ' s\n')
    f.close()
    tracemalloc.stop()
