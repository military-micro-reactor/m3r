from math import pi, floor, sin, cos, hypot
from sys import stdout
from random import random
import numpy as np

from particle import Particle

compact_r = 0.5  # cm
compact_h = 1.2  # cm
_pf = 0.40

default_particle = Particle()


def random_fill(pv, cr, ch, pf):
    """
    Fill compact with particles by random placement.
    Terminates after 50K attempts.
    Will not fully fill over 0.3 packing fraction due to naive algorithm.
    Parameters:
    pv: avg. particle volume
    cr: compact radius
    ch: compact height
    pf: packing fraction of particles
    """
    cv = pi * (compact_r**2) * compact_h
    n = floor(pf * cv / pv)
    print("Filling compact (r={}, h={}) with {} particles".format(cr, ch, n))
    print()

    # I know it's bad.
    # Try to randomly fill the compact :(
    # Eventually, a transition to a repeated structures model
    #    will guarantee better performance and take less space
    particle_list = []
    i = 0
    while len(particle_list) < n and i < 100000:
        stdout.write("\r")
        next_p = Particle()
        try_r = (cr - next_p.total_radius) * random()
        try_z = next_p.total_radius + (compact_h -
                                       2 * next_p.total_radius) * random()
        try_theta = 2 * pi * random()
        next_p.set_origin(
            cos(try_theta) * try_r,
            sin(try_theta) * try_r, try_z)
        intersects = False
        for p in particle_list:
            if p.intersects(next_p):
                intersects = True
                break
        i += 1
        if not intersects:
            stdout.write("Placing particle {:10} of {}".format(
                len(particle_list) + 1, n))
            particle_list.append(next_p)
    print()
    print("{} particles placed before termination.".format(len(particle_list)))
    print("Actual packing fraction is {:.4f}.".format(
        (pv * len(particle_list)) / cv))
    return particle_list


def _cell_intersects(compact_rad, cell_x, cell_y, cell_rad):
    for i in [-1, 1]:
        for j in [-1, 1]:
            if compact_rad < hypot(cell_x + i * cell_rad,
                                   cell_y + j * cell_rad):
                return True
    return False


def generate_fill_lattice(pr, pv, compact_rad, compact_height, pf):
    """
    Fills a compact with a lattice for URAN card particle modeling
    Paramters:
    pr: particle radius
    pv: particle volume
    compact_rad: compact radius
    compact_height: compact height
    pf: packing fraction of particles. Must be less than roughly .40
    Returns:
    The number of cells in one plane dimensions,
        the number of cells in the axial dimension,
        and a grid of ones and zeros representing where there are cells that can hold particles
    """
    print(
        "Filling compact (r = {}, h = {}) with particles (r = {:.4f}) at packing fraction {:.4f}"
        .format(compact_rad, compact_height, pr, pf))
    cv = pi * (compact_r**2) * compact_h
    num_particles = floor(pf * cv / pv)
    print("Attempting to place {} particles".format(num_particles))
    print()
    cell_rad = 0.5 * compact_rad  # 5% "wiggle room" for the URAN card
    tries = 0
    while tries < 1000:

        cell_dia = 2 * cell_rad
        # number of cells that can fit into one compact radius
        #   minus the radius of the center cell
        num_radial_cells = floor((compact_rad - cell_rad) / (cell_dia))
        num_plane_cells = 2 * num_radial_cells + 1
        num_axial_cells = floor(compact_height / cell_dia)
        cell_pass_grid = np.ones((num_plane_cells, num_plane_cells))
        # pass grid is misaligned with (x, y) coordinate space
        # could use meshgrid, but a simple subtraction is easier
        start_xy = -1 * (num_radial_cells) * cell_dia

        for i in range(num_plane_cells):
            for j in range(num_plane_cells):
                if _cell_intersects(compact_rad, start_xy + i * cell_dia,
                                    start_xy + j * cell_dia, cell_rad):
                    cell_pass_grid[i][j] = 0

        total_cells = np.sum(cell_pass_grid) * num_axial_cells
        stdout.write("Mesh cells: {} \r".format(total_cells))

        if total_cells > num_particles:
            print(cell_rad)
            return num_plane_cells, num_axial_cells, cell_pass_grid

        # there has to be a better way to do this
        cell_rad = 0.99 * cell_rad
        if cell_rad < pr:
            print(
                "\nUnable to fill with cubic lattice at this packing fraction")
            print("actual packing fraction: {:.4f}".format(
                (num_particles * pv) / cv))
            return 0, None, None
        tries += 1
    return 0, None, None


if __name__ == '__main__':
    random_fill(default_particle.volume(), compact_r, compact_h, 0.25)
    print(
        generate_fill_lattice(default_particle.total_radius,
                              default_particle.volume(), compact_r, compact_h,
                              0.30))
