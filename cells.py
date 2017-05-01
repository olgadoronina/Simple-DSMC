from params import *
from velocity import *


class Cells(object):
    def __init__(self):
        self.all = [[Subcell(i, j) for j in range(n_subdiv ** 2)] for i in range(N_x * N_y)]

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
                        N_in_subcell[i * n_subdiv + l, j * n_subdiv + m] = len(
                            self.all[j * N_x + i][m * n_subdiv + l].set)
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
            return int(round(N_coll))
        else:
            return 0

    def collision_pairs(self):
        """Return set of pairs chosen for collision"""
        pairs = set()
        counter = 0
        N_coll = self.calc_N_coll()
        # print('N_coll = ', N_coll)
        for i in range(N_coll):
            mol1, mol2 = rand.sample(self.set, 2)
            c_r = sqrt((mol1.u-mol2.u)**2+(mol1.v-mol2.v)**2+(mol1.w-mol2.w)**2)
            self.c_r_max = max(c_r, self.c_r_max)
            if c_r / self.c_r_max > rand.random():
                pairs.add((mol1, mol2, c_r))
                counter+=1
        # print('N_actual', counter)
        return pairs

    def collision(self):
        """Perform collisions and assign new velocities"""

        def c_r_new_direction(c_r):
            cos_theta = 2*rand.random()-1
            sin_theta = sqrt(1 - cos_theta ** 2)
            phi = rand.uniform(0, 2 * pi)
            return np.array([c_r*cos_theta, c_r*sin_theta * cos(phi), c_r*sin_theta * sin(phi)])

        pairs = self.collision_pairs()
        # print(len(pairs))
        if len(pairs) == 0:
            return
        for (mol1, mol2, c_r) in self.collision_pairs():
            # p=0
            # if norm(mol1.velocity())>1500 or norm(mol2.velocity())>1500:
            #     p = 1
            #     print('c_r old', c_r)
            #     print('old', norm(mol1.velocity()), norm(mol2.velocity()),norm(mol1.velocity())+norm(mol2.velocity()))
            #     print(norm((mol1.velocity() + mol2.velocity()) / 2))
            v_m = (mol1.velocity() + mol2.velocity()) / 2.
            c_r_new = c_r_new_direction(c_r)
            mol1.new_velocity(v_m - 0.5 * c_r_new)
            mol2.new_velocity(v_m + 0.5 * c_r_new)
            # if p:
            #     print('new', norm(mol1.velocity()), norm(mol2.velocity()),norm(mol1.velocity())+norm(mol2.velocity()) )
            #     print('new', norm(mol1.velocity() - mol2.velocity()))

        # print('c_r_max', self.c_r_max)