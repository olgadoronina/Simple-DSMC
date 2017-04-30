from params import *
from velocity import *

class Cells(object):
    def __init__(self):
        self.all = [[Subcell(i,j) for j in range(n_subdiv ** 2)] for i in range(N_x * N_y)]

    def add_to_subcell(self, molecule):
        """Add the molecules to the cells and subcells"""
        index_cell, index_subcell = molecule.get_cell()
        self.all[index_cell][index_subcell].set.add(molecule)

    def calcNumber(self):
        N_in_cell = np.zeros((N_x, N_y), dtype=float)
        N_in_subcell = np.zeros((N_x * n_subdiv, N_y * n_subdiv), dtype=float)
        for i in range(N_x):
            for j in range(N_y):
                for l in range(n_subdiv):
                    for m in range(n_subdiv):
                        N_in_cell[i, j] += len(self.all[j * N_x + i][m * n_subdiv + l].set)
                        N_in_subcell[i * n_subdiv + l, j * n_subdiv + m] = len(self.all[j * N_x + i][m * n_subdiv + l].set)
        print('mean N in cells = ', np.mean(N_in_cell))
        print('mean N in subcells = ', np.mean(N_in_subcell))
        return (N_in_cell, N_in_subcell)

    def calcTemperature(self):
        T_in_cell = np.zeros((N_x, N_y), dtype=float)
        T_in_subcell = np.zeros((N_x * n_subdiv, N_y * n_subdiv), dtype=float)
        for i in range(N_x):
            for j in range(N_y):
                for l in range(n_subdiv):
                    for m in range(n_subdiv):
                        for molecule in self.all[j * N_x + i][m * n_subdiv + l].set:
                            T_in_cell[i, j] += molecule.temperature()
                            T_in_subcell[i * n_subdiv + l, j * n_subdiv + m] += molecule.temperature()

        return (T_in_cell, T_in_subcell)


class Subcell(object):
    def __init__(self, index_cell, index_subcell):
        self.set = set()
        self.index_subcell = index_subcell
        self.index_cell = index_cell
        self.c_r_max = c_r_max_initial

    def calc_N_coll(self):
        """Calculate N of molecules in subcell to collide"""
        N_in = len(self.set)
        # print('N_in', N_in)
        if N_in > 1:
            N_coll = 1 / (2 * V_subcell) * N_in * (N_in - 1) * Fn * sigma * self.c_r_max * time_step
            # print("N_coll = ", round(N_coll), N_in)
            return min(round(N_coll), (factorial(N_in) / factorial(2) / factorial(N_in-2)))
        else:
            return 0

    def collision_pairs(self):
        """Return set of pairs chosen for collision"""
        pairs = set()
        counter = 0
        N_coll = self.calc_N_coll()
        while counter < N_coll:
            mol1, mol2 = rand.sample(self.set,2)
            c_r = mol1.speed() - mol2.speed()
            self.c_r_max = max(c_r, self.c_r_max)
            # print('c_r_max', self.c_r_max)
            if c_r / self.c_r_max > rand.random():
                counter += 1
                pairs.add((mol1, mol2, c_r))
        return pairs

    def collision(self):
        """Perform collisions and assign new velocities"""

        def c_r_new_direction():
            theta = rand.uniform(0,2*pi)
            phi  = rand.uniform(0,2*pi)
            return np.array([cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi)])

        for (mol1, mol2, c_r) in self.collision_pairs():
            v_m = (mol1.velocity()+ mol2.velocity())/2
            c_r_new = c_r*c_r_new_direction()
            mol1.new_velocity(v_m + 0.5 * c_r_new)
            mol2.new_velocity(v_m - 0.5 * c_r_new)
