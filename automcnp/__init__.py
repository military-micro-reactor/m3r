# Initialization for automcnp module

# testing
from .particle import Particle
k = Particle()
print(k.mcnp_kernel_material(1000))
print(k.mcnp_transform(1000))
print(k.mcnp_surfaces(100))
print(k.mcnp_cells(100, 100))
