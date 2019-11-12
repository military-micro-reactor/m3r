from math import pi, floor, sin, cos
from random import random
from automcnp.particle import Particle
from automcnp.defaults import NOMINAL_PACKING_FRAC as _npf


particle_list = []

compact_r = 1.0  # cm
compact_h = 1.2  # cm
compact_vol = pi*(compact_r**2)*compact_h

default_particle = Particle()

n_particles = floor(compact_vol/default_particle.volume())

print("Filling compact (r={}, h={}) with {} particles".format(compact_r,
                                                              compact_h,
                                                              n_particles))
i = 0
while len(particle_list) < n_particles and i<100000:
    next_p = Particle()
    try_r = (compact_r - next_p.total_radius)*random()
    try_z = next_p.total_radius + (compact_h - 2*next_p.total_radius)*random()
    try_theta = 2*pi*random()
    next_p.set_origin(cos(try_theta)*try_r, sin(try_theta)*try_r, try_z)
    intersects = False
    for p in particle_list:
        if p.intersects(next_p):
            print("Intersection. Skipping this particle.")
            intersects = True
            break
    i += 1
    if not intersects:
        particle_list.append(next_p)
print(len(particle_list))
    
