from math import pi as _PI
from math import hypot

from defaults import *


# MEAN KERNEL PARAMETERS -- AGR-2 KERNELS
_u_enrichment = 0.20  # uranium enrichment wt%
_kernel_density = 11.0  # fuel kernel density g/cc
_kernel_fraction_o = 0.0861  # wt%
_kernel_fraction_c = 0.0178  # wt%
_kernel_fraction_238u = 0.7169  # wt%
_kernel_fraction_235u = 0.1792  # wt%
# KERNEL PARAMETER FROM AGR5+ SPEC
_kernel_radius = 425 * MICRON


class Particle:
    """
    Models a TRISO particle as a series of concentric spheres.
    Defines methods to to auto-generate MCNP cells, surfaces, and materials.
    """

    def __init__(self,
                 kernel_rad=_kernel_radius,
                 kernel_density=_kernel_density,
                 kernel_spec=((8016, _kernel_fraction_o),
                              (6012, _kernel_fraction_c),
                              (92238, _kernel_fraction_238u),
                              (92235, _kernel_fraction_235u)),
                 buffer_density=BUFFER_DENSITY,
                 buffer_thickness=BUFFER_THICKNESS,
                 ipyc_density=IPYC_DENSITY,
                 ipyc_thickness=IPYC_THICKNESS,
                 sic_density=SIC_DENSITY,
                 sic_thickness=SIC_THICKNESS,
                 opyc_density=OPYC_DENSITY,
                 opyc_thickness=OPYC_THICKNESS,
                 use_weight_frac=True):
        self.kernel_rad = kernel_rad
        self.kernel_spec = kernel_spec
        self.kernel_density = kernel_density
        self.buffer_density = buffer_density
        self.buffer_thickness = buffer_thickness
        self.ipyc_density = ipyc_density
        self.ipyc_thickness = ipyc_thickness
        self.sic_density = sic_density
        self.sic_thickness = sic_thickness
        self.opyc_density = opyc_density
        self.opyc_thickness = opyc_thickness

        self.total_radius = (
            self.kernel_rad + self.buffer_thickness + self.ipyc_thickness +
            self.sic_thickness + self.opyc_thickness)
        # dx, dy, dz used to build MCNP cell transform
        self.dx = 0
        self.dy = 0
        self.dz = 0

    def set_origin(self, x, y, z):
        """
        Sets the origin of the particle to the ordered triple (x, y, z)
        Paramters:
        x, y, z
        """
        self.dx = x
        self.dy = y
        self.dz = z

    def volume(self):
        return _PI * (4 / 3) * (self.total_radius)**3

    def intersects(self, other):
        """
        Check to see if this particle intersects another particle.
        Parameters:
        other: a particle object
        Returns: 
        True if the two particles intersect.
        """
        return (hypot(
            hypot(self.dx - other.dx, self.dy - other.dy), self.dz - other.dz)
                <= (self.total_radius + other.total_radius))

    def mcnp_transform(self, transform_n):
        """
        Generates the MCNP transformation card for this particle.
        Parameters:
        transform_n: The transformation number

        Returns:
        The MCNP transformation card in cosine format
        """
        # TRN dx dy dz 1 0 0 0 1 0 0 0 1 1
        transform_card = "TR{N} {dx} {dy} {dz} 1 0 0 0 1 0 0 0 1 1"
        self.transform_n = transform_n
        return transform_card.format(
            N=transform_n, dx=self.dx, dy=self.dy, dz=self.dz)

    def mcnp_kernel_material(self, material_n):
        self.kernel_material = material_n
        material_card = "{prefix}{zaid} {density}"
        mcard = material_card.format(
            prefix="m{} ".format(self.kernel_material),
            zaid=self.kernel_spec[0][0],
            density=MCNP_DENSITY_FACTOR * self.kernel_spec[0][1])
        for i in range(1, len(self.kernel_spec)):
            mcard = mcard + material_card.format(
                prefix="\n      ",
                zaid=self.kernel_spec[i][0],
                density=MCNP_DENSITY_FACTOR * self.kernel_spec[i][1])
        return mcard

    def mcnp_surfaces(self, surface_n, transform_n=0):
        """
        Generates a list of MCNP surface cards for this particle.
        Each particle requires 5 surfaces.
        Paramters:
        surface_n: the surface number to start numbering the particle surfaces
        (optional) transform_n: the particle transformation number. 
        Will use self.transform_n if the mcnp_transform method has been called.
        """
        if transform_n == 0:
            try:
                transform_n = self.transform_n
            except AttributeError:
                print("TRANSFORMATION NOT YET DEFINED")
                print("USING NO TRANSFORMATION")
                transform_n = ""
        self.surfaces = [
            surface_n, surface_n + 1, surface_n + 2, surface_n + 3,
            surface_n + 4
        ]
        surf_card = "{n} {tr} so {r:0.8f}"
        surfaces = []
        surfaces.append(
            surf_card.format(
                n=self.surfaces[0], tr=transform_n, r=self.kernel_rad))
        surfaces.append(
            surf_card.format(
                n=self.surfaces[1],
                tr=transform_n,
                r=(self.kernel_rad + self.buffer_thickness)))
        surfaces.append(
            surf_card.format(
                n=self.surfaces[2],
                tr=transform_n,
                r=(self.kernel_rad + self.buffer_thickness +
                   self.ipyc_thickness)))
        surfaces.append(
            surf_card.format(
                n=self.surfaces[3],
                tr=transform_n,
                r=(self.kernel_rad + self.buffer_thickness +
                   self.ipyc_thickness + self.sic_thickness)))
        surfaces.append(
            surf_card.format(
                n=self.surfaces[4],
                tr=transform_n,
                r=(self.kernel_rad + self.buffer_thickness +
                   self.ipyc_thickness + self.sic_thickness +
                   self.opyc_thickness)))
        return surfaces

    def mcnp_cells(self,
                   cell_n,
                   universe_n,
                   buffer_mat_n=100,
                   ipyc_mat_n=101,
                   sic_mat_n=102,
                   opyc_mat_n=103):

        self.cells = [cell_n, cell_n + 1, cell_n + 2, cell_n + 3]
        cell_card = "{n} {m} {rho} {geom} imp:n=1 u={u}"
        cells = []
        cells.append(
            cell_card.format(
                n=self.cells[0],
                m=self.kernel_material,
                rho=MCNP_DENSITY_FACTOR * self.kernel_density,
                geom=-1 * self.surfaces[0],
                u=universe_n))
        materials = [buffer_mat_n, ipyc_mat_n, sic_mat_n,
                     opyc_mat_n]
        densities = [self.buffer_density, self.ipyc_density,
                   self.sic_density, self.opyc_density]
        for i in range(1, 4):
            cells.append(
                cell_card.format(
                    n=self.cells[i],
                    m=materials[i-1],
                    rho=MCNP_DENSITY_FACTOR * densities[i-1],
                    geom="{} {}".format(-1 * self.surfaces[i],
                                        self.surfaces[i - 1]),
                    u=universe_n))
        return cells
