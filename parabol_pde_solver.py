# CODE FOR ADVANCE TOPICS, SOLVES A PARABOLIC PARTIAL DIFFERENTIAL EQUATION ON A 2D SQUARE DOMAIN

# IMPORT OF LIBRARIES
import numpy as np  # arrays and numerical application
from scipy import sparse  # sparse matrices
from scipy.sparse import linalg

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

    topology = topology.astype(int)

    return nodes, topology

def const_voronoi_cells(top, area, nds):
    """Calculates the voronoi cells associated to each element assuming triangular regular mesh with
    all elements being of area area"""
    vor_cells = np.zeros(len(nds))
    for i in range(len(nds)):
        vor_cells[i] = np.count_nonzero(top == i+1)*area/3

    return vor_cells


def hf_calc(nds, top, area, vor_cells):
    """Function used to assembly the P and H matrices element wise using the elements in
    top and the nodes in nds. The parameter sv is the vector containing the S values for each element"""

    # prepare rows and col values for the sparse matrix
    rows = top.flatten()
    rows = rows - np.ones(len(rows))
    rows = rows.astype(int)
    col = (np.repeat(np.arange(1, len(top) + 1, 1), 3))
    col = col - np.ones(len(col))
    col = col.astype(int)
    val = np.repeat(np.nan,len(rows))

    temp = sparse.csr_matrix((val, (col, rows)))
    # memory pre allocation
    H = np.transpose(temp) * temp
    P = np.transpose(temp) * temp
    fvt = np.zeros(len(nds))
    # This will be a point to split, how many processes can I use? How should I split?
    for nel, elem in enumerate(topology):

        h_loc = local_H(nds, elem, area)
        f_loc = local_f(elem, vor_cells, nds)
        for i in range(3):
            ind_i = elem[i] - 1

            fvt[ind_i] = fvt[ind_i] + f_loc[i]

            for j in range(3):
                ind_j = elem[j] - 1
                if np.isnan(H[ind_i,ind_j]):
                    H[ind_i, ind_j] =  h_loc[i, j]
                else:
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

    loc_f = np.zeros(shape=(3,1))

    for i,ni in enumerate([ni, nj, nk]):
        loc_f[i] = f_fun(nds[ni-1,:]) * vor_cells[ni-1]
    return loc_f


def dirichlet_bc(H, f, uv, dbc_nodes, penalty = 10**20):
    """Function to add the Dirichlet BCs to the problem: H is the sparse stiffness matrix, f the vector of known values,
    xv and yv the x and y coordinates of the nodes respectively and u the values of the solution"""

    for node in dir_nodes:
        H[node, node] = penalty
        f[node] = uv[node]*penalty

    return f, H


def u_fun(node):
    """Function to calculate the true value of the unknown function on the node node needed for the boundary
    conditions """
    x = node[0]
    y = node[1]

    return  x**2 + y**2 - (x**2)*(y**2) - 1 #np.cos(2 * x ** 2 - y) + y ** 2


def f_fun(node):
    """Function to calculate the true value of the unknown function on the node node needed for the boundary
    conditions """
    x = node[0]
    y = node[1]

    return 4 - 2*x**2 - 2*y**2 #-(-4 * np.sin(2 * x ** 2 - y) - (16 * x ** 2 + 1) * np.cos(2 * x ** 2 - y) + 2)


# MAIN CODE
extr = 10
step = 0.05


x = np.arange(-extr, extr+step, step)
y = np.arange(-extr, extr+step, step)

side = len(x)

nodes, topology = create_simple_mesh(x, y)

shape = np.array([[1, nodes[topology[0, 0] - 1, 0], nodes[topology[0, 0] - 1, 1]],
                  [1, nodes[topology[0, 1] - 1, 0], nodes[topology[0, 1] - 1, 1]],
                  [1, nodes[topology[0, 2] - 1, 0], nodes[topology[0, 2] - 1, 1]]])

area = np.linalg.det(shape) / 2
vor_cells = const_voronoi_cells(topology, area, nodes)

uv_true = np.apply_along_axis(u_fun, 1, nodes)

H, fv = hf_calc(nodes, topology, area, vor_cells)

dir_nodes = np.array(np.concatenate((range(side), range( side * (side-1), side**2, 1), range(side, side*(side-1), side), range(2*side-1, side*(side-1), side))))

fv, H = dirichlet_bc(H, fv, uv_true, dir_nodes, penalty=10**20)
u0 = np.zeros(len(uv_true))
uv, info = linalg.cg(H, fv, atol=10**(-5), tol=10**(-5))





