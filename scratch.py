# CODE FOR ADVANCE TOPICS, SOLVES A PARABOLIC PARTIAL DIFFERENTIAL EQUATION ON A 2D SQUARE DOMAIN

# IMPORT OF LIBRARIES
import numpy as np  # arrays and numerical application
from scipy import sparse  # sparse matrices


# FUNCTIONS


def create_simple_mesh(x, y, s=1):
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

    sv = np.repeat(s, len(topology))
    topology = topology.astype(int)

    return nodes, topology, sv


def phf_calc(nds, top, area, sv, fv):
    """Function used to assembly the P and H matrices element wise using the elements in
    top and the nodes in nds. The parameter sv is the vector containing the S values for each element"""

    # prepare rows and col values for the sparse matrix
    rows = top.flatten()
    rows = rows - np.ones(len(rows))
    rows = rows.astype(int)
    col = (np.repeat(np.arange(1, len(top) + 1, 1), 3))
    col = col - np.ones(len(col))
    col = col.astype(int)
    val = np.zeros(shape=(len(rows)))

    temp = sparse.csr_matrix((val, (col, rows)))
    # memory pre allocation
    H = np.transpose(temp) * temp
    P = np.transpose(temp) * temp
    fvt = np.zeros(len(fv))
    # This will be a point to split, how many processes can I use? How should I split?
    for nel, elem in enumerate(topology):

        h_loc = local_H(nds, elem, area)
        p_loc = local_P(nel, elem, area, sv)
        f_loc = local_f(fv, elem, area)
        for i in range(3):
            ind_i = elem[i] - 1

            fvt[ind_i] = fvt[ind_i] + f_loc[i]

            for j in range(3):
                ind_j = elem[j] - 1

                H[ind_i, ind_j] = H[ind_i, ind_j] + h_loc[i, j]
                P[ind_i, ind_j] = P[ind_i, ind_j] + p_loc[i, j]

    return H, P, fvt


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


def local_P(el_ind, elem, area, sv):
    """Function used to calculate the local capacity matrix on an element"""

    ni = elem[0]
    nj = elem[1]
    nk = elem[2]

    loc_p = np.zeros(shape=(3, 3))

    for i in range(3):
        for j in range(3):
            if i == j:
                loc_p[i, j] = sv[el_ind]*area/6
            else:
                loc_p[i, j] = sv[el_ind]*area/12

    return loc_p


def local_f(f, elem, area):
    """Function to calculate the actual value of f element-wise"""

    ni = elem[0]
    nj = elem[1]
    nk = elem[2]

    loc_f = np.zeros(shape=(3,1))

    for i,ni in enumerate([ni, nj, nk]):
        loc_f[i] = f[ni-1] * (area / 3)

    return loc_f


def dirichlet_bc(H, f, xv, yv, uv):
    """Function to add the Dirichlet BCs to the problem: H is the sparse stiffness matrix, f the vector of known values,
    xv and yv the x and y coordinates of the nodes respectively and u the values of the solution"""

    return ()


def neumann_bc(f, xv, yv, qv, nds):
    """Function to add the Neumann BCs to the problem: f is the vector of known values, xv and yv are the x and y
    coordinates of the nodes respectively, qv are the values of the q function on the boundaries of the domain and nds
    are the nodes of the mesh"""
    return ()


def pcg_jacobi(H, f, u0, epsi=10 ** (-10), kmax=1000):
    """Function that implements the preconditioned conjugate gradient method for solving linear systems,
    H is the stiffness matrix of the problem, f the vector of known values, epsilon the tolerance of the error,
    kmax the maximum number of iterations allowed and u0 the vector containing the starting point"""
    return ()


def theta_method(H, P, u_time, f, theta, t_step, x_val, y_val, uv, qv, nds, epsi=10 ** (-10), kmax=1000):
    """Function used to calculate the finite difference solution of the partial differential equation using the theta
    method. H is the stiffness matrix of the problem, P the capacity matrix, u_time the 2D array containing all the
    computed solutions (needed for each subsequent step), f the vector of known values, theta the value of theta for the
    theta method and t_step the time step interval to use. The remaining parameter (xv, yv, uv, qv, nds, epsilon and
    kmax) are the parameters needed for the pcg_jacobi and neumann_bc functions"""
    return ()


def u_fun(node):
    """Function to calculate the true value of the unknown function on the node node needed for the boundary
    conditions """
    x = node[0]
    y = node[1]

    return np.cos(2 * x ** 2 - y) + y ** 2


def f_fun(node):
    """Function to calculate the true value of the unknown function on the node node needed for the boundary
    conditions """
    x = node[0]
    y = node[1]

    return -(-4 * np.sin(2 * x ** 2 - y) - (16 * x ** 2 + 1) * np.cos(2 * x ** 2 - y) + 2)


def q_fun(node):
    """Function to calculate the true value of the unknown function on the node node needed for the boundary
    conditions """
    x = node[0]
    y = node[1]

    return -np.sin(2 * x ** 2 - y) * 4 * x


# MAIN CODE


# TESTING FUNCTIONS

x = np.arange(-0.5, 0.6, 1)
y = np.arange(-0.5, 0.6, 1)

nodes, topology, sv = create_simple_mesh(x, y, 10)

shape = np.array([[1, nodes[topology[0, 0] - 1, 0], nodes[topology[0, 0] - 1, 1]],
                  [1, nodes[topology[0, 1] - 1, 0], nodes[topology[0, 1] - 1, 1]],
                  [1, nodes[topology[0, 2] - 1, 0], nodes[topology[0, 2] - 1, 1]]])

area = np.linalg.det(shape) / 2

uv_true = np.apply_along_axis(u_fun, 1, nodes)
fv = np.apply_along_axis(f_fun, 1, nodes)

H, P, fv = phf_calc(nodes, topology, area, sv, fv)

dir_nodes = [[range(len(x))]]