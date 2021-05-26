# CODE FOR ADVANCE TOPICS, SOLVES A PARABOLIC PARTIAL DIFFERENTIAL EQUATION ON A 2D SQUARE DOMAIN

# IMPORT OF LIBRARIES
import numpy as np


# FUNCTIONS
def create_simple_mesh(x, y, s):
    """Function used to create a simple mesh of a squared domain"""
    side = len(x)
    xv = np.tile(x, side)  # X coordinates in row by row
    yv = np.repeat(y, side)  # Y coordinates in column*row
    nodes = np.asarray([xv, yv]).transpose()

    topology = np.zeros(shape=(((side-1)**2)*2, 3))
    ind = 0
    nd = 1

    while nd <= len(nodes):  # cycle over nodes
        if nd % side != 0 and nd <= side * (side - 1):  # skip nodes on the right edge
            topology[ind, :] = [nd, nd + side + 1, nd + side] # top of upper triang
            topology[ind + 1, :] = [nd, nd + 1, nd + side + 1] # top of lower triang
            ind = ind + 2 # two rows at a time
        nd = nd + 1 # go to next node

    sv = np.repeat(s, len(topology))

    return nodes, topology, sv


def ph_calc(nds, top, sv):
    """Function used to assembly the P and H matrices element wise using the elements in
    top and the nodes in nds. The parameter sv is the vector containing the S values for each element"""
    return ()


def f_calc(f, nds, top):
    """Function to calculate the actual value of f to be used in the computations,
    f is the initial f (analytical solution), nds the list of nodes an top the topology of the mesh"""
    return ()


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


# MAIN CODE


# TESTING FUNCTIONS

x = np.arange(-1, 1.5, 0.5)
y = np.arange(-1, 1.5, 0.5)
print(y)
print(x)

nodes, topology, sv = create_simple_mesh(x, y, 5)
print(nodes)
print(topology)
print(sv)
