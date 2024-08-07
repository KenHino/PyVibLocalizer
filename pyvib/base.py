from collections import defaultdict, Counter
import copy
import itertools
import math
import mendeleev
import numpy as np

from typing import List, Dict, Optional, Tuple, Union

from pyvib.local_cls import Localizer, GroupLocalizer
import pyvib.units as units

def dist(A: List[float], B: List[float]):
    vecA = np.array(A)
    vecB = np.array(B)
    return np.linalg.norm(vecA - vecB)

def read_minfo(file_name: str, use_trans: Optional[bool] =False,
        use_rot: Optional[bool] =False
        ) -> Tuple[List[Union[str, Tuple[float, float, float]]],
                List[float], List[List[float]]]:
    geom = []
    freq = []
    disp = []
    read_geom = False
    read_disp = False
    read_freq = False
    file = open(file_name,"r")
    for i, line in enumerate(file):
        if line == '[ Atomic Data ]\n':
            read_geom = True
            continue

        elif line == 'Translational Frequency\n' and use_trans:
            read_freq = True
            read_disp = False
            continue

        elif line == 'Translational vector\n' and use_trans:
            read_freq = False
            read_disp = True
            continue

        elif line == 'Rotational Frequency\n' and use_rot:
            read_freq = True
            read_disp = False
            continue

        elif line == 'Rotational vector\n' and use_rot:
            read_freq = False
            read_disp = True
            continue

        elif line == 'Vibrational Frequency\n':
            read_freq = True
            read_disp = False
            continue

        elif line == 'Vibrational vector\n':
            read_freq = False
            read_disp = True
            continue

        elif read_geom:
            if line == '\n':
                read_geom = False
                continue
            else:
                line = line.replace(',', '').replace('E', 'e')
                words = line.split()
                if len(words) == 6:
                    geom.append([words[0], [float(w) for w in words[3:6]]])
                else:
                    natom = int(words[0])
                continue

        elif read_freq:
            if line[-2] == ' ':
                line = line.replace(',', '').replace('e', 'e')
                words = line.split()
                freq.extend([float(w) for w in words])
                continue

        elif read_disp:
            if line[-2] == ' ':
                line = line.replace(',', '').replace('e', 'e')
                words = line.split()
                disp.extend([float(w) for w in words])
                continue

    disp = [disp[3*natom*k:3*natom*(k+1)] for k in range(len(disp)//3//natom)]

    return (geom, freq, disp)

def write_minfo(file_name: str, freq: List[float], disp: List[List[float]]):
    nmode = len(freq)
    contents = 'Vibrational Frequency\n'
    contents += str(nmode) + '\n'
    for imode in range(nmode):
        contents += "{:>15.8e}".format(freq[imode])
        if imode % 5 == 4 or imode == nmode -1:
            contents += ' \n'
        else:
            contents += ', '

    contents += 'Vibrational vector\n'
    for imode in range(nmode):
        contents += 'Mode ' + str(imode+1) + '\n'
        contents += str(len(disp[imode])) + '\n'
        for k in range(len(disp[imode])):
            contents += "{:>15.8e}".format(disp[imode][k])
            if k % 5 == 4 or k == len(disp[imode]) - 1:
                contents += ' \n'
            else:
                contents += ', '

    file = open(file_name, 'w')
    file.write(contents)
    file.close()

def read_fchk_g16(file_path : str
    ) -> Tuple[str, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read Hessian in `.fchk` files from Gaussian16.

    Args:
        file_path (str): Path to ``fchk`` file.
        use_trans (bool): Use translational vector. Defaults to ``False``.
        use_rot (bool): Use rotatinal vector. Defaults to ``False``.

    Returns:
        Tuple: multiple contents listed below.

    Return contents:
        - List[List[Union[str, Tuple[float, float, float]]]]: Geometry
        - numpy.ndarray: mass-weighted hessian.

    """
    from mendeleev import element
    atom = []
    d1_hess_cartesian = []
    mass = []
    geom = []

    atom_flag = False
    geom_flag = False
    hessian_flag = False
    mass_flag = False

    file = open(file_path,"r")
    lines = file.readlines()
    file.close()
    for i, line in enumerate(lines):
        if line.startswith('Cartesian Force Constants'):
            hessian_flag=True
            continue
        elif line.startswith('Vib-AtMass'):
            mass_flag = True
            continue
        elif line.startswith('Opt point       1 Geometries'):
            geom_flag = True
            continue
        elif line.startswith('Atomic numbers'):
            atom_flag = True
            continue

        if hessian_flag | mass_flag | geom_flag | atom_flag:
            if line[-5] != 'E' and not atom_flag:
                hessian_flag = False
                mass_flag = False
                geom_flag = False
                continue
            if line.startswith('Nuclear charges'):
                atom_flag = False
                continue

            if hessian_flag:
                d1_hess_cartesian.extend(list(map(float, line.split())))
            elif mass_flag:
                mass.extend(list(map(float, line.split())))
            elif geom_flag:
                geom.extend(list(map(float, line.split())))
            elif atom_flag:
                atomic_num = list(map(int, line.split()))
                atom.extend([element(n).symbol for n in atomic_num])


    def to_symmetric(v,size):
        out = np.zeros((size, size), v.dtype)
        out[np.tri(size, dtype=bool)] = v
        return out.T + out - np.diag(np.diag(out))

    natom = len(atom)
    coord= np.array(geom).reshape(natom, 3)
    hess_cartesian = to_symmetric(np.array(d1_hess_cartesian), natom*3)
    mass = np.array(mass)
    trans_mass_weighted = 1 / np.sqrt(np.repeat(mass * units.DALTON, 3)).reshape(-1,1)
    mass_weighted_hessian = np.multiply(hess_cartesian, trans_mass_weighted@trans_mass_weighted.T)

    if False:
        E, V = scipy.linalg.eigh(mass_weighted_hessian)
        freq = np.sqrt(E) * const.au_in_cm1

        ndof = V.shape[0]
        disp_mwc = np.array(V.T).reshape(ndof, natom, 3)
        disp = np.zeros_like(disp_mwc)
        for iatom in range(natom):
            disp[:,iatom, :] = disp_mwc[:,iatom,:] / np.sqrt(mass[iatom])
    geom = [[a,tuple(c)] for a,c in zip(atom, coord)]
    return (geom, mass_weighted_hessian)

class Vibration:
    r""" Main Vibraion class

    This class allows localization and visualization of vibrational modes.

    Attributes:
       unit (str): Unit of displacement vector. Defaults to ``Bohr``.
       freq (Optional, List[float]): The list of frequencies in cm-1.
       mw_disp (Optional, List[List[float]]): The list of mass-weighted displacement vectors.
                `disp[imode][iatom_xyz]` gives displacement from refernce coordinates.
       mw_hess (Optional, np.ndarray): mass-weighted hessian
       geom (List[List[str, Tuple[float,flotat,float]]]): Referenece geometry.
            ``geom[iatom][0]`` gives element symbol.
            ``geom[iatom][1]`` gives the refernce cartesian coordinate.
       bond (List[List[Tuple[int,int]]]): The pair of atom coordinates which bonded each other.
       nmode (int): Number of modes.
       natom (int): Number of atoms.
       Q_mat (numpy.ndarray): Unitary matrix aligned displacement vectors.

    Args:
       geom (List[List[str, Tuple[float,float,float]]]): Referenece geometry.
            ``geom[iatom][0]`` gives element symbol.
            ``geom[iatom][1]`` gives the refernce cartesian coordinate.
       freq (List[float]): The list of frequencies in cm-1.
       disp (List[List[float]]): The list of displacement vectors.
                ``disp[imode][iatom_xyz]`` gives displacement from refernce coordinates.

    """
    def __init__(self,
        geom: List[List[Union[str, Tuple[float,float,float]]]],
        freq: Optional[List[float]] =None,
        mw_disp: Optional[List[List[float]]] =None,
        mw_hess: Optional[np.ndarray] =None,
        unit_xyz: Optional[str] ='Bohr',
        unit_xyz_hess: Optional[str] ='Bohr',
        unit_omega: Optional[str] ='Hartree',
        unit_mass: Optional[str] ='a.u.',
        disp: Optional[List[List[float]]] =None):

        self.atom = [a[0] for a in geom]
        
        if disp is not None:
            self.disp = disp
            mw_disp = self.set_mwdisp()
        elif mw_disp is None:
            mw_disp = [[0.0 for _ in range(len(geom)*3)]]
            freq = [-1.0]
        else:
            assert len(freq) == len(mw_disp)
        

        if unit_xyz.lower() in ['bohr', 'au', 'a.u.']:
            pass
        elif unit_xyz.lower() in ['angstrom']:
            geom = [[g[0], (np.array(g[1]) * units.ANGSTROM).tolist()] for g in geom]
            mw_disp = (np.array(mw_disp) * units.ANGSTROM).tolist()
        else:
            raise NotImplementedError

        if unit_omega.lower() in ['cm1', 'cm', 'cm-1', 'kayser']:
            freq = (np.array(freq) * units.CM1).tolist()
        elif unit_omega.lower() in ['ev']:
            freq = (np.array(freq) * units.EV).tolist()
        elif unit_omega.lower() in ['au', 'a.u.', 'hartree']:
            pass
        else:
            raise NotImplementedError

        if unit_mass.lower() in ['amu', 'dalton']:
            mw_disp = (np.array(mw_disp) * math.sqrt(units.DALTON)).tolist()
        elif unit_mass.lower() in ['a.u.', 'au', 'emu']:
            pass
        else:
            raise NotImplementedError

        self.freq = np.array(freq)
        self.mw_disp = np.array(mw_disp)

        self.geom = geom
        self.coord = [c[1] for c in geom]
        self.bond = [[self.coord[i], self.coord[k]] for (i,k) in itertools.combinations(range(len(self.atom)), 2) \
                if (mendeleev.element(self.atom[i]).covalent_radius_pyykko +
                mendeleev.element(self.atom[k]).covalent_radius_pyykko)*1.0e-02 + 0.1>
                dist(self.coord[i], self.coord[k]) / units.ANGSTROM]

        self.natom = len(self.geom)
        self.set_disp()
        self.Q_mat = np.matrix(self.disp).T

        if mw_hess is not None:
            self.diagonalize(mw_hess=mw_hess,
            unit_xyz=unit_xyz_hess, unit_mass=unit_mass, unit_omega=unit_omega)

        self.nmode = len(self.mw_disp)

    def set_disp(self):
        self.disp = np.array(self.mw_disp)
        for i, a in enumerate(self.atom):
            self.disp[:, 3*i:3*i+3] /= math.sqrt(mendeleev.element(a).atomic_weight * units.DALTON)
        return self.disp
    
    def set_mwdisp(self):
        self.mw_disp = np.array(self.disp)
        for i, a in enumerate(self.atom):
            self.mw_disp[:, 3*i:3*i+3] *= math.sqrt(mendeleev.element(a).atomic_weight * units.DALTON)
        return self.mw_disp

    def visualize(self, arrow_scale: Optional[float]=100,
            blender: Optional[bool] =False, atom_number: Optional[bool]=False):
        """Visualize vibrational modes.

        Vibration mode visualization with Matplotlib and tinker.

        Args:
            arrow_scale (float): The scale of displacement arrows. Defaults to ``1``.
            blender (Optional[bool]): View by blender.
            atom_number (Optional[bool]): Plot atom number

        Examples:
            >>> sim = Vibration(geom, freq, disp)
            >>> sim.visualize()

        """
        self.blender = blender
        self.atom_number = atom_number
        if self.blender:
            from pyvib.blender import Visualizer
            vis = Visualizer(self, arrow_scale)
            vis.show()

        else:
            from mpl_toolkits.mplot3d import Axes3D
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            import matplotlib.pyplot as plt
            import tkinter as tk
            from tkinter import ttk
            from pyvib.visual_cls import Visualizer
            ### root object ###
            root = tk.Tk()
            root.title("PyVibVisualizer")

            inputFrame = ttk.Frame(root)
            buttonFrame = ttk.Frame(root)
            graphFrame = ttk.Frame(root)

            ### inputFrame ###
            inputData = Visualizer(inputFrame, self, arrow_scale)
            inputFrame.pack()

            ### buttonFrame ###
            ButtonWidth = 10
            UpdateButton = tk.Button(buttonFrame, text="View", width=ButtonWidth, \
                command = lambda : inputData.plot(Canvas, ax, fig))
            UpdateButton.grid(row = 0, column = 0)
            buttonFrame.pack()

            #graph initialize
            fig = plt.figure(figsize=(8,8))
            ax = fig.add_subplot(111, projection='3d')
            ax.set_box_aspect((1,1,1))
            ax.axis('off')
            Canvas = FigureCanvasTkAgg(fig, master = graphFrame)
            Canvas.get_tk_widget().pack()
            graphFrame.pack()

            #continue
            root.mainloop()

    def localize(self, option: Optional[str] ='Pipek-Mezey',
            window :Optional[float]=400):
        r"""Localize vibrational modes.

        Vibration mode localization with local_cls module.

        - Pipek-Mezey metric

        $$
            \xi_{\mathrm{at}}\left(\widetilde{\boldsymbol{Q}}^{\mathrm{sub}}\right)=\sum_{p=1}^{k} \sum_{i=1}^{n}\left(\tilde{C}_{i p}^{\mathrm{sub}}\right)^{2}
            \\
            \tilde{C}_{i p}^{\mathrm{sub}}=\sum_{\alpha=x, y, z}\left(\tilde{Q}_{i \alpha, p}^{\mathrm{sub}}\right)^{2}
        $$

        - Boys metric

        $$
            \xi_{\mathrm{dist}}\left(\widetilde{\boldsymbol{Q}}^{\mathrm{sub}}\right)=\sum_{p=1}^{k}\left(\boldsymbol{R}_{p}^{\text {center }}\right)^{2}
            \\
            \boldsymbol{R}_{p}^{\text {center }}=\sum_{i=1}^{n} \tilde{C}_{i p}^{\mathrm{sub}} \boldsymbol{R}_{i}
        $$

        Args:
            option (str): The metric of localization. ``Boys`` or ``Pipek-Mezey``.
            window (float): window frequency in cm-1.

        Examples:
            >>> sim = Vibration(geom, freq, disp)
            >>> sim.localize()

        """
        loc = Localizer(self, option, window)
        self.Q_mat, self.freq = loc.run()
        for x in range(self.nmode):
            for y in range(self.natom*3):
                self.disp[x][y] = self.Q_mat[y,x]

        print(option,'metric = ',loc.metric(self.Q_mat))
        return (self.disp, self.freq)

    def diagonalize(self,
        mw_hess: np.ndarray,
        unit_xyz: Optional[str] ='au',
        unit_omega: Optional[str] ='au',
        unit_mass: Optional[str] ='au'
        ) -> Tuple[List[List[float]], List[float]]:
        mw_disp, freq = self.group_localize(mw_hess, domain=[[k for k in range(self.natom)]], 
                                            unit_xyz=unit_xyz, unit_omega=unit_omega, unit_mass=unit_mass)
        self.freq = freq[6:]
        self.mw_disp = mw_disp[6:]
        self.disp = self.disp[6:]
        return (self.disp, self.freq)

    def group_localize(self,
        mw_hess: np.ndarray,
        domain: List[List[float]],
        unit_xyz: Optional[str] ='au',
        unit_omega: Optional[str] ='au',
        unit_mass: Optional[str] ='au'
        ) -> Tuple[List[List[float]], List[float]]:
        r"""Generate Group Localized Coordinate from mass-weighted hessian

        Args:
            mwhess (np.ndarray): mass-weighted hessian
            domain (List[List[float]]): atomic domain such as [[0,1,2],[3,4]]

        """
        mw_hess = np.array(mw_hess)
        if unit_omega.lower() in ['cm1', 'cm', 'cm-1', 'kayser']:
            '''
            E = 1/2 m omega^2 x^2
            hess = omega^2
            mwc = sqrt{m}x
            '''
            mw_hess *= units.CM1**2
        elif unit_omega.lower() in ['ev']:
            mw_hess *= units.EV**2
        elif unit_omega.lower() in ['au', 'a.u.', 'hartree']:
            pass
        else:
            raise NotImplementedError

        loc = GroupLocalizer(self, mw_hess, domain)
        self.Q_mat, freq = loc.run()
        self.freq = [0.0 if math.isnan(f) else f for f in freq]            

        if unit_xyz.lower() in ['bohr', 'au', 'a.u.']:
            pass
        elif unit_xyz.lower() in ['angstrom']:
            self.Q_mat *= units.ANGSTROM
        else:
            raise NotImplementedError

        if unit_mass.lower() in ['amu', 'dalton']:
            self.Q_mat *= math.sqrt(units.DALTON)
        elif unit_mass.lower() in ['a.u.', 'au', 'emu']:
            pass
        else:
            raise NotImplementedError
        self.mw_disp = self.Q_mat.T.tolist()
        self.set_disp()
        return (self.disp, self.freq)
