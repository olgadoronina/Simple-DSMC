from params import *
from cells import *
from velocity import *


class Molecule(object):
    def __init__(self, x=None, y=None, flux=None):

        def acceptance_rejection_uniform(beta):
            """Acceptance-rejection procedure for gaussian distribution in equilibrium gas"""
            PP_max = lambda u: np.exp(-beta ** 2 * u ** 2)
            a, b = -3 / beta, 3 / beta  # cut-off for pdf
            while True:
                R_f = rand.random()
                u = rand.uniform(a, b)
                if R_f < PP_max(u):
                    return u

        def init_coords(x, y):
            """Calculate initial conditions for sample molecule
            using uniform distribution"""
            if x == None:
                self.x = rand.uniform(0, X_max)
            else:
                self.x = x
            if y == None:
                self.y = rand.uniform(0, Y_max)
            else:
                self.y = y
            if self.x == 0:
                self.x += TINY
            if self.y == 0:
                self.y += TINY

        def init_velocity(beta):
            """Calculate initial velocities
            using Acceptance-Rejection method for u'
            and direct pair method for v' and w' """
            phi = rand.uniform(0, 2 * pi)
            rho = sqrt(-log(rand.random())) / beta
            if not flux:
                self.u = acceptance_rejection_uniform(beta)
            elif flux == 'input':
                self.u = acceptance_rejection(s, beta)
            elif flux == 'output':
                self.u = -acceptance_rejection(-s, beta)
            else:
                sys.exit("Molecule(): init_velocity(): flux undefined")

            # self.v = rho * cos(phi)
            # self.w = rho * sin(phi)
            self.v = acceptance_rejection_uniform(beta)
            self.w = acceptance_rejection_uniform(beta)

        init_coords(x, y)
        init_velocity(beta)

    def get_cell(self):
        """Return index of cell and subcell where molecule is.
        If statements are needed to count molecules, which on the border"""
        if self.x == X_max:
            L = N_x - 1
            i = n_subdiv - 1
        else:
            L = int((self.x) / cell_size)
            i = int((self.x - L * cell_size) / subcell_size)
        if self.y == Y_max:
            M = N_y - 1
            j = n_subdiv - 1
        else:
            M = int((self.y) / cell_size)
            j = int((self.y - M * cell_size) / subcell_size)
        index_cell = M * N_x + L
        index_subcell = j * n_subdiv + i
        return index_cell, index_subcell

    def update_coord(self):
        """Update the coords of molecule after time step"""

        def new_coord(time_left=time_step):
            """Calculate the new coordinates after time step
            without wall consideration"""
            x_new = self.x + self.u * time_left
            y_new = self.y + self.v * time_left
            intersection_flag = 0
            if x_new < 0 or x_new > X_max or y_new < 0 or y_new > Y_max:
                intersection_flag = 1
            else:
                self.x, self.y = x_new, y_new
            return intersection_flag

        def intersection(time_left):
            """Define the coordinates of intersection,
            return the time left after intersection"""
            t = time_left
            bound = None
            for x_bound in ['left', 'right']:
                time = (boundary[x_bound] - self.x) / self.u
                if time > 0 and time <= t:
                    t = time
                    bound = x_bound
            for y_bound in ['down', 'up']:
                time = (boundary[y_bound] - self.y) / self.v
                if time > 0 and time <= t:
                    t = time
                    bound = y_bound
            if bound is 'left':
                self.x = boundary[bound] + TINY
                self.y += self.v * t
            if bound is 'right':
                self.x = boundary[bound] - TINY
                self.y += self.v * t
            if bound is 'down':
                self.x += self.u * t
                self.y = boundary[bound] + TINY
            if bound is 'up':
                self.x += self.u * t
                self.y = boundary[bound] - TINY
            if bound == None:
                print(self.x, self.y, self.u, self.v, time, time_left)
                sys.exit("Molecule:update_coord(self):intersection(): bound is undefined")
            return time_left - t, bound

        flag_bound_collision = new_coord()
        time_left = time_step
        while flag_bound_collision:
            time_left, bound_of_intersection = intersection(time_left)
            if INPUT and bound_of_intersection in ['left', 'right']:
                self.x = 'delete'
                break
            self.u, self.v, self.w = newVelocity(bound_of_intersection, [self.u, self.v, self.w])
            if abs(time_left) < TINY:
                break
            flag_bound_collision = new_coord(time_left)

    def speed(self):
        return sqrt(self.u ** 2 + self.v ** 2 + self.w ** 2)

    def velocity(self):
        return np.array([self.u, self.v, self.w])

    def new_velocity(self, vel):
        self.u = vel[0]
        self.v = vel[1]
        self.w = vel[2]

    def temperature(self):
        speed_sqr = self.u ** 2 + self.v ** 2 + self.w ** 2
        return m * speed_sqr / 3 / k
