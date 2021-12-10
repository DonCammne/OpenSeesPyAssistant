# Module with the fibers
#   Carmine Schipani, 2021

from numpy.matrixlib import mat
from openseespy.opensees import *
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Wedge
import openseespy.postprocessing.ops_vis as opsv
import numpy as np
import os
import math
from copy import copy, deepcopy
from OpenSeesPyHelper.Section import *
from OpenSeesPyHelper.DataManagement import *
from OpenSeesPyHelper.ErrorHandling import *
from OpenSeesPyHelper.Units import *
from OpenSeesPyHelper.MaterialModels import *

class Fibers(DataManagement):
    pass

class FibersRect(Fibers):
    # Class that stores funcions and material properties of the rectangular fiber section.
    # Warning: the units should be m and N
    # Coordinates: plotting coordinte (x, y) = fiber section coordinate (-z, y)
    
    def __init__(self, ID: int, b, d, Ay, D_hoops, e, unconf_mat_ID: int, conf_mat_ID: int, bars_mat_ID: int,
        bars_x: np.ndarray, ranges_y: np.ndarray, discr_core: list, discr_cover_lateral: list, discr_cover_topbottom: list, GJ = 0.0):

        # Check
        if ID < 1: raise NegativeValue()
        if b < 0: raise NegativeValue()
        if d < 0: raise NegativeValue()
        if Ay < 0: raise NegativeValue()
        if D_hoops < 0: raise NegativeValue()
        if e < 0: raise NegativeValue()
        if unconf_mat_ID < 1: raise NegativeValue()
        if conf_mat_ID < 1: raise NegativeValue()
        if bars_mat_ID < 1: raise NegativeValue()
        if np.size(bars_x) != np.size(ranges_y)-1: raise WrongDimension()
        geometry_tol = 5*mm_unit
        for bars in bars_x:
            if abs(np.sum(bars) - b) > geometry_tol: raise InconsistentGeometry()
        if abs(np.sum(ranges_y) - d) > geometry_tol: raise InconsistentGeometry()
        if e > b/2 or e > d/2: raise InconsistentGeometry()
        if len(discr_core) != 2: raise WrongDimension()
        if len(discr_cover_lateral) != 2: raise WrongDimension()
        if len(discr_cover_topbottom) != 2: raise WrongDimension()
        if GJ < 0: raise NegativeValue()

        # Arguments
        self.ID = ID
        self.b = b
        self.d = d
        self.Ay = Ay
        self.D_hoops = D_hoops
        self.e = e
        self.unconf_mat_ID = unconf_mat_ID
        self.conf_mat_ID = conf_mat_ID
        self.bars_mat_ID = bars_mat_ID
        self.bars_x = deepcopy(bars_x)
        self.ranges_y = copy(ranges_y)
        self.discr_core = discr_core
        self.discr_cover_lateral = discr_cover_lateral
        self.discr_cover_topbottom = discr_cover_topbottom
        self.GJ = GJ #TODO: check if with 0, no problems or how to compute it

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        self.ReInit()

    def ReInit(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Memebers
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        
        # Parameters
        z1 = self.b/2
        y1 = self.d/2
        zc = z1-self.e-self.D_hoops
        yc = y1-self.e-self.D_hoops

        # Create the concrete core fibers
        core = [-yc, -zc, yc, zc]
        core_cmd = ['patch', 'rect', self.conf_mat_ID, self.discr_core[1], self.discr_core[0], *core]

        # Create the concrete cover fibers (top, bottom, left, right)
        cover_up = [yc, -z1, y1, z1]
        cover_down = [-y1, -z1, -yc, z1]
        cover_left = [-yc, zc, yc, z1]
        cover_right = [-yc, -z1, yc, -zc]
        cover_up_cmd = ['patch', 'rect', self.unconf_mat_ID, self.discr_cover_topbottom[1], self.discr_cover_topbottom[0], *cover_up]
        cover_down_cmd = ['patch', 'rect', self.unconf_mat_ID, self.discr_cover_topbottom[1], self.discr_cover_topbottom[0], *cover_down]
        cover_left_cmd = ['patch', 'rect', self.unconf_mat_ID, self.discr_cover_lateral[1], self.discr_cover_lateral[0], *cover_left]
        cover_right_cmd = ['patch', 'rect', self.unconf_mat_ID, self.discr_cover_lateral[1], self.discr_cover_lateral[0], *cover_right]
        self.fib_sec = [core_cmd, cover_up_cmd, cover_down_cmd, cover_left_cmd, cover_right_cmd]
        
        # Create the reinforcing fibers (top, middle, bottom)
        nr_bars = 0
        for range in self.bars_x:
            nr_bars += np.size(range)-1
        rebarY = -np.cumsum(self.ranges_y[0:-1]) + y1
        self.rebarYZ = np.zeros((nr_bars, 2))

        iter = 0
        for ii, Y in enumerate(rebarY):
            rebarZ = -np.cumsum(self.bars_x[ii][0:-1]) + z1
            for Z in rebarZ:
                self.rebarYZ[iter, :] = [Y, Z]
                iter = iter + 1
        
        for YZ in self.rebarYZ:
            self.fib_sec.append(['layer', 'bar', self.bars_mat_ID, self.Ay, *YZ])

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "FibersRect"], # Tag for differentiating different data
            ["ID", self.ID],
            ["section_name_tag", self.section_name_tag],
            ["b", self.b],
            ["d", self.d],
            ["Ay", self.Ay],
            ["D_hoops", self.D_hoops],
            ["e", self.e],
            ["GJ", self.GJ],
            ["conf_mat_ID", self.conf_mat_ID],
            ["discr_core", self.discr_core],
            ["unconf_mat_ID", self.unconf_mat_ID],
            ["discr_cover_topbottom", self.discr_cover_topbottom],
            ["discr_cover_lateral", self.discr_cover_lateral],
            ["bars_mat_ID", self.bars_mat_ID],
            ["bars_x", self.bars_x],
            ["ranges_y", self.ranges_y]]

    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the material model (optional).
        """
        print("")
        print("Requested info for FibersRect, ID = {}".format(self.ID))
        print("Section associated: {} ".format(self.section_name_tag))
        print("Base b = {} mm and depth d = {} mm".format(self.b/mm_unit, self.d/mm_unit))
        print("Confined material model ID = {}".format(self.conf_mat_ID))
        print("Unconfined material model ID = {}".format(self.unconf_mat_ID))
        print("Bars material model ID = {}".format(self.bars_mat_ID))
        print("Discretisation in the core [x dir, y dir] = {}".format(self.discr_core))
        print("Discretisation in the lateral covers [x dir, y dir] = {}".format(self.discr_cover_lateral))
        print("Discretisation in the top and bottom covers [x dir, y dir] = {}".format(self.discr_cover_topbottom))
        print("")

        if plot:
            plot_fiber_section(self.fib_sec, ['#808080', '#D3D3D3', 'k'])
            
            if block:
                plt.show()


    def CreateFibers(self):
        create_fiber_section(self.fib_sec)

class FibersRectRCRectShape(FibersRect):
    def __init__(self, ID: int, ele: RCRectShape, unconf_mat_ID: int, conf_mat_ID: int, bars_mat_ID: int,
        discr_core: list, discr_cover_lateral: list, discr_cover_topbottom: list, GJ=0):
        super().__init__(ID, ele.b, ele.d, ele.Ay, ele.D_hoops, ele.e, unconf_mat_ID, conf_mat_ID, bars_mat_ID,
            ele.bars_position_x, ele.bars_ranges_position_y, discr_core, discr_cover_lateral, discr_cover_topbottom, GJ=GJ)
        self.section_name_tag = ele.name_tag
        self.UpdateStoredData()


def plot_fiber_section(fiber_info, fill_shapes = True, matcolor=['#808080', '#D3D3D3', 'r', 'b', 'g', 'y']):
    """Plot fiber cross-section.

    Args:
        fiber_info (list): list of lists in the format similar to the parameters
            for the section, layer (circ, straight or bar), patch, fiber OpenSees commands (coordinate sys = yzlocal!!!)

        fill_shapes (bool): True - filled fibers with color specified in matcolor
            list, False - no color, only the outline of fibers

        matcolor (list): sequence of colors for various material tags
            assigned to fibers (number of colors should not be less than the number of mat model used)

    Examples:
        ::

            fib_sec_1 = [['section', 'Fiber', 1, '-GJ', 1.0e6],
                         ['patch', 'quad', 1, 4, 1,  0.032, 0.317, -0.311, 0.067, -0.266, 0.005, 0.077, 0.254],  # noqa: E501
                         ['patch', 'quad', 1, 1, 4,  -0.075, 0.144, -0.114, 0.116, 0.075, -0.144, 0.114, -0.116],  # noqa: E501
                         ['patch', 'quad', 1, 4, 1,  0.266, -0.005,  -0.077, -0.254,  -0.032, -0.317,  0.311, -0.067]  # noqa: E501
                         ]
            opsv.fiber_info_to_cmds(fib_sec_1)
            matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
            opsv.plot_fiber_section(fib_sec_1, matcolor=matcolor)
            plt.axis('equal')
            # plt.savefig('fibsec_rc.png')
            plt.show()

    Notes:
        ``fiber_info`` can be reused by means of a python helper function
            ``ops_vis.fiber_info_to_cmds(fiber_info_1)``

    See also:
        ``ops_vis.fiber_info_to_cmds()``
    """

    mat_to_col = {}
    fig, ax = plt.subplots()
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.grid(False)

    for item in fiber_info:

        if item[0] == 'layer':
            matID = item[2]
            mat_to_col = __assignColorToMat(matID, mat_to_col, matcolor)
            if item[1] == 'bar':
                As = item[3]
                Iy = item[4]
                Iz = item[5]
                r = np.sqrt(As / np.pi)
                bar = Circle((-Iz, Iy), r, ec='k', fc='k', zorder=10)
                ax.add_patch(bar)
            if item[1] == 'straight':
                n_bars = item[3]
                As = item[4]
                Iy, Iz, Jy, Jz = item[5], item[6], item[7], item[8]
                r = np.sqrt(As / np.pi)
                Y = np.linspace(Iy, Jy, n_bars)
                Z = np.linspace(Iz, Jz, n_bars)
                for zi, yi in zip(Z, Y):
                    bar = Circle((-zi, yi), r, ec='k', fc=mat_to_col[matID], zorder=10)
                    ax.add_patch(bar)
            if item[1] == 'circ':
                n_bars, As = item[3], item[4]
                yC, zC, arc_radius = item[5], item[6], item[7]
                if len(item) > 8:
                    a0_deg, a1_deg = item[8], item[9]
                else:
                    a0_deg, a1_deg = 0., 360. - 360./n_bars

                a0_rad, a1_rad = np.pi * a0_deg / 180., np.pi * a1_deg / 180.
                r_bar = np.sqrt(As / np.pi)
                thetas = np.linspace(a0_rad, a1_rad, n_bars)
                Y = yC + arc_radius * np.cos(thetas)
                Z = zC + arc_radius * np.sin(thetas)
                for zi, yi in zip(Z, Y):
                    bar = Circle((-zi, yi), r_bar, ec='k', fc=mat_to_col[matID], zorder=10)
                    ax.add_patch(bar)

        if (item[0] == 'patch' and (item[1] == 'quad' or item[1] == 'quadr' or
                                  item[1] == 'rect')):
            matID, nIJ, nJK = item[2], item[3], item[4]
            mat_to_col = __assignColorToMat(matID, mat_to_col, matcolor)


            if item[1] == 'quad' or item[1] == 'quadr':
                Iy, Iz, Jy, Jz = item[5], item[6], item[7], item[8]
                Ky, Kz, Ly, Lz = item[9], item[10], item[11], item[12]

            if item[1] == 'rect':
                Iy, Iz, Ky, Kz = item[5], item[6], item[7], item[8]
                Jy, Jz, Ly, Lz = Ky, Iz, Iy, Kz

            # check for convexity (vector products)
            outIJxIK = (Jy-Iy)*(Kz-Iz) - (Ky-Iy)*(Jz-Iz)
            outIKxIL = (Ky-Iy)*(Lz-Iz) - (Ly-Iy)*(Kz-Iz)
            # check if I, J, L points are colinear
            outIJxIL = (Jy-Iy)*(Lz-Iz) - (Ly-Iy)*(Jz-Iz)
            # outJKxJL = (Ky-Jy)*(Lz-Jz) - (Ly-Jy)*(Kz-Jz)

            if outIJxIK <= 0 or outIKxIL <= 0 or outIJxIL <= 0:
                print('Warning! Patch quad is non-convex or counter-clockwise defined or has at least 3 colinear points in line')

            IJz, IJy = np.linspace(Iz, Jz, nIJ+1), np.linspace(Iy, Jy, nIJ+1)
            JKz, JKy = np.linspace(Jz, Kz, nJK+1), np.linspace(Jy, Ky, nJK+1)
            LKz, LKy = np.linspace(Lz, Kz, nIJ+1), np.linspace(Ly, Ky, nIJ+1)
            ILz, ILy = np.linspace(Iz, Lz, nJK+1), np.linspace(Iy, Ly, nJK+1)

            if fill_shapes:
                Z = np.zeros((nIJ+1, nJK+1))
                Y = np.zeros((nIJ+1, nJK+1))

                for j in range(nIJ+1):
                    Z[j, :] = np.linspace(IJz[j], LKz[j], nJK+1)
                    Y[j, :] = np.linspace(IJy[j], LKy[j], nJK+1)

                for j in range(nIJ):
                    for k in range(nJK):
                        zy = np.array([[-Z[j, k], Y[j, k]],
                                       [-Z[j, k+1], Y[j, k+1]],
                                       [-Z[j+1, k+1], Y[j+1, k+1]],
                                       [-Z[j+1, k], Y[j+1, k]]])
                        poly = Polygon(zy, True, ec='k', fc=mat_to_col[matID])
                        ax.add_patch(poly)

            else:
                # horizontal lines
                for az, bz, ay, by in zip(IJz, LKz, IJy, LKy):
                    plt.plot([-az, -bz], [ay, by], 'b-', zorder=1)

                # vertical lines
                for az, bz, ay, by in zip(JKz, ILz, JKy, ILy):
                    plt.plot([-az, -bz], [ay, by], 'b-', zorder=1)

        if item[0] == 'patch' and item[1] == 'circ':
            matID, nc, nr = item[2], item[3], item[4]
            mat_to_col = __assignColorToMat(matID, mat_to_col, matcolor)

            yC, zC, ri, re = item[5], item[6], item[7], item[8]
            a0, a1 = item[9], item[10]

            dr = (re - ri) / nr
            dth = (a1 - a0) / nc

            for j in range(nr):
                rj = ri + j * dr
                rj1 = rj + dr

                for i in range(nc):
                    thi = a0 + i * dth
                    thi1 = thi + dth
                    wedge = Wedge((yC, -zC), rj1, thi, thi1, width=dr, ec='k',
                                  lw=1, fc=mat_to_col[matID])
                    ax.add_patch(wedge)

    ax.axis('equal')


def __assignColorToMat(matID: int, mat_to_col: dict, matcolor: list):
    if not matID in mat_to_col:
        if len(mat_to_col) >= len(matcolor):
            print("Warning: not enough colors defined for fiber section plot (white used)")
            mat_to_col[matID] = 'w'
        else:
            mat_to_col[matID] = matcolor[len(mat_to_col)]
    return mat_to_col


def create_fiber_section(fiber_info):
    """Reuses fiber_info to define fiber section in OpenSees.

    At present it is not possible to extract fiber section data from
    the OpenSees domain, this function is a workaround. The idea is to
    prepare data similar to the one the regular OpenSees commands
    (``section('Fiber', ...)``, ``fiber()``, ``patch()`` and/or
    ``layer()``) require.

    Args:
        fiber_info (list): is a list of fiber section data. First sub-list
        also defines the torsional stiffness (GJ).

    Warning:

    If you use this function, do not issue the regular OpenSees:
    section, Fiber, Patch or Layer commands.

    See also:

    ``ops_vis.plot_fiber_section()``

    """
    for dat in fiber_info:
        if dat[0] == 'section':
            secTag, GJ = dat[2], dat[4]
            section('Fiber', secTag, '-GJ', GJ)

        if dat[0] == 'layer':
            matTag = dat[2]
            if dat[1] == 'straight':
                n_bars = dat[3]
                As = dat[4]
                Iy, Iz, Jy, Jz = dat[5], dat[6], dat[7], dat[8]
                layer('straight', matTag, n_bars, As, Iy, Iz, Jy, Jz)
            if dat[1] == 'bar':
                As = dat[3]
                Iy = dat[4]
                Iz = dat[5]
                layer('straight', matTag, 1, As, Iy, Iz, Iy, Iz)

        if dat[0] == 'patch':
            matTag = dat[2]
            nIJ = dat[3]
            nJK = dat[4]

            if dat[1] == 'quad' or dat[1] == 'quadr':
                Iy, Iz, Jy, Jz = dat[5], dat[6], dat[7], dat[8]
                Ky, Kz, Ly, Lz = dat[9], dat[10], dat[11], dat[12]
                patch('quad', matTag, nIJ, nJK, Iy, Iz, Jy, Jz, Ky, Kz,
                        Ly, Lz)

            if dat[1] == 'rect':
                Iy, Iz, Ky, Kz = dat[5], dat[6], dat[7], dat[8]
                Jy, Jz, Ly, Lz = Ky, Iz, Iy, Kz
                patch('rect', matTag, nIJ, nJK, Iy, Iz, Ky, Kz)

        
