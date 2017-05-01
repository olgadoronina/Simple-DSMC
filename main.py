#!/usr/bin/python3
# -*- coding: utf-8 -*-
from plots import *
from molecules import *
from params import *
from cells import *
import multiprocessing as mp

## Global structure
cells = Cells()

def moving_molecules(cell):
    """Calculate new coordinates for all molecules after one time step"""
    for subcell in cell:
        for molecule in subcell.set:
            molecule.update_coord()
    return cell


def put_in_subcells(cell):
    """Redistribute molecules in subcells according new coordinates"""
    for subcell in cell:
        for molecule in list(subcell.set):
            if molecule.x == 'delete':
                subcell.set.remove(molecule)
            else:
                index_cell, index_subcell = molecule.get_cell()
                if index_cell != subcell.index_cell or index_subcell != subcell.index_subcell:
                    subcell.set.remove(molecule)
                    cells.all[index_cell][index_subcell].set.add(molecule)
    return cell

def collision(cell):
    """perform molecular collision in each subcell"""
    for subcell in cell:
        subcell.collision()
    return cell

def main():
    print('lambda/3 = ', lambd / 3, '\tsubcell size = ', subcell_size)
    print('N molecules per subcell = ', Ne / N_x / N_y / n_subdiv ** 2)
    print('n initially ', n_initial)
    print('sigma = ', sigma)
    print('time step = ', time_step)
    print('N real = ', N, '\tF_N = ', Fn)
    print('Ne = ', Ne)
    print('V = ', V)
    print('N input', N_input)
    print('N output', N_output)
    if PARALLEL:
        print("\n%d workers" % mp.cpu_count())

    start = time()
    #######################################################################
    ## Initial conditions
    #
    # if INPUT:
    #     for i in range(N_input):
    #         molecule = Molecule(x=0, y=rand.uniform(0, Y_max), flux='input')
    #         cells.add_to_subcell(molecule)
    #     for i in range(N_output):
    #         molecule = Molecule(x=X_max, y=rand.uniform(0, Y_max), flux='output')
    #         cells.add_to_subcell(molecule)
    # else:
    for i in range(Ne):
        molecule = Molecule()
        cells.add_to_subcell(molecule)
    N_in_cell, N_in_subcell = cells.calcNumber()
    plotDistribution(N_in_subcell, np.linspace(np.min(N_in_subcell), np.max(N_in_subcell), 21), 'Subcells', "nearest", 0)
    plotDistribution(N_in_cell, np.linspace(np.min(N_in_cell), np.max(N_in_cell), 21), 'Cells', "nearest", 0)
    #######################################################################


    #######################################################################
    ## Main Loop
    for step in range(N_steps):
        start_step = time()
        print("Time step #", step)
        # ########################################################################
        ## Moving molecules
        start_motion = time()
        if PARALLEL:
            p = mp.Pool(mp.cpu_count())
            all_cells = p.map(moving_molecules, cells.all, 1)  # molecule.update_coord() in parallel
            p.close()
            p.join()
            cells.all = all_cells
        else:
            for cell in cells.all:
                moving_molecules(cell)
        end_motion = time()
        timer(start_motion, end_motion, 'Time of motion')
        # ########################################################################
        ## Put molecules in cells
        start_cell = time()
        for cell in cells.all:
            put_in_subcells(cell)
        end_cell = time()
        timer(start_cell, end_cell, 'Time of putting molecules in cells')
        # ########################################################################
        ## collision
        if COLLISION:
            start_coll = time()
            if PARALLEL:
                p = mp.Pool(mp.cpu_count())
                all_cells = p.map(collision, cells.all)  # molecule.update_coord() in parallel
                p.close()
                p.join()
                cells.all = all_cells
            else:
                for cell in cells.all:
                    collision(cell)
            end_coll = time()
            timer(start_coll, end_coll, 'Time of molecular collision ')
        N_in_cell, N_in_subcell = cells.calcNumber()
        plotDistribution(N_in_subcell, np.linspace(np.min(N_in_subcell), np.max(N_in_subcell), 21), 'Subcells', "nearest", 0)
        plotDistribution(N_in_cell, np.linspace(np.min(N_in_cell), np.max(N_in_cell), 21), 'Cells', "nearest", 1)

        T_in_cell, T_in_subcell = cells.calcTemperature()
        T_in_cell = np.divide(T_in_cell, N_in_cell, out=np.zeros_like(T_in_cell), where=N_in_cell != 0)
        T_in_subcell = np.divide(T_in_subcell, N_in_subcell, out=np.zeros_like(T_in_subcell), where=N_in_subcell != 0)
        print('mean T in cells = ', np.mean(T_in_cell))
        print('mean T in subcells = ', np.mean(T_in_subcell))
        plotDistribution(T_in_subcell, np.linspace(np.min(T_in_subcell), np.max(T_in_subcell), 21), 'Subcells', "nearest", 0)
        plotDistribution(T_in_cell, np.linspace(np.min(T_in_cell), np.max(T_in_cell), 21), 'Cells', "nearest", 1)
        # ########################################################################
        ## input and output
        if INPUT:
            start_input = time()
            for i in range(N_input):
                molecule = Molecule(x=0, y=rand.uniform(0, Y_max), flux='input')
                cells.add_to_subcell(molecule)
            for i in range(N_output):
                molecule = Molecule(x=X_max, y=rand.uniform(0, Y_max), flux='output')
                cells.add_to_subcell(molecule)
            end_input = time()
            timer(start_input, end_input, 'Time per input')
        end_step = time()
        timer(start_step, end_step, 'Time per whole step')

    # ########################################################################


    # ########################################################################
    N_in_cell, N_in_subcell = cells.calcNumber()
    plotDistribution(N_in_subcell, np.linspace(np.min(N_in_subcell), np.max(N_in_subcell), 21), 'Subcells', "nearest", 0)
    plotDistribution(N_in_cell, np.linspace(np.min(N_in_cell), np.max(N_in_cell), 21), 'Cells', "nearest", 1)

    T_in_cell, T_in_subcell = cells.calcTemperature()
    T_in_cell = np.divide(T_in_cell, N_in_cell, out=np.zeros_like(T_in_cell), where=N_in_cell != 0)
    T_in_subcell = np.divide(T_in_subcell, N_in_subcell, out=np.zeros_like(T_in_subcell), where=N_in_subcell != 0)
    print('mean T in cells = ', np.mean(T_in_cell))
    print('mean T in subcells = ', np.mean(T_in_subcell))
    plotDistribution(T_in_subcell, np.linspace(np.min(T_in_subcell), np.max(T_in_subcell), 21), 'Subcells', "nearest", 0)
    plotDistribution(T_in_cell, np.linspace(np.min(T_in_cell), np.max(T_in_cell), 21), 'Cells', "nearest", 1)

    end = time()
    timer(start, end, 'Time for whole simulation')

    plot_T_function(N_in_cell)
    plot_N_function(T_in_cell)

    print(N_in_cell[0, 0], N_in_subcell[0, 0])
    print(N_in_cell[0, 50], N_in_subcell[0, 100])
    print(N_in_cell[50, 0], N_in_subcell[100, 0])
    print(N_in_cell[50, 50], N_in_subcell[100, 100])
    print(T_in_cell[0, 0], T_in_subcell[0, 0])
    print(T_in_cell[0, 50], T_in_subcell[0, 100])
    print(T_in_cell[50, 0], T_in_subcell[100, 0])
    print(T_in_cell[50, 50], T_in_subcell[100, 100])


if __name__ == "__main__":
    main()
