import itertools
from collections import defaultdict, Counter
import copy
import math
try:
    import tkinter as tk
    from tkinter import ttk
except:
    import tk
    from tk import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from typing import List, Dict, Optional, Tuple, Union
import mendeleev

from visual_cls import Visualizer
from local_cls import Localizer, GroupLocalizer
import units


def dist(A, B):
    dum = 0.0
    for a, b in zip(A, B):
        dum += (a-b)**2
    return math.sqrt(dum)


def read_minfo(file_name, use_trans=False, use_rot=False):
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

    return geom, freq, disp


def write_minfo(file_name, freq, disp):
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


class Vibration:
    """ Main Vibraion class

    This class allows localization and visualization of vibrational modes.

    Attributes:
        unit (str) : Unit of displacement vector. Defaults to ``Bohr``.
        freq (List(float)) : The list of frequencies in cm-1.
        disp (List[List[float]]) : The list of displacement vectors. \
                ``disp[imode][iatom_xyz]`` gives displacement from refernce coordinates.
        geom (List[List[str, Tuple[float,flotat,float]]]) : Referenece geometry. \
            ``geom[iatom][0]`` gives element symbol. \
            ``geom[iatom][1]`` gives the refernce cartesian coordinate.
        bond (List[List[Tuple[int,int]]]]) : The pair of atom coordinates which bonded each other.
        nmode (int) : Number of modes.
        natom (int) : Number of atoms.
        Q_mat (numpy.ndarray) : Unitary matrix aligned displacement vectors.

    Args:
        geom (List[List[str, Tuple[float,float,float]]]) : Referenece geometry. \
            ``geom[iatom][0]`` gives element symbol. \
            ``geom[iatom][1]`` gives the refernce cartesian coordinate.
        freq (List[float]) : The list of frequencies in cm-1.
        disp (List[List[float]]) : The list of displacement vectors. \
                ``disp[imode][iatom_xyz]`` gives displacement from refernce coordinates.
    """
    def __init__(self, 
        geom : List[List[Union[str, Tuple[float,float,float]]]], 
        freq : List[float] =None, 
        disp_mw : List[List[float]] =None, 
        unit_xyz : Optional[str] ='Bohr', 
        unit_freq : Optional[str] ='cm1', 
        unit_mass : Optional[str] ='AMU'):
        if disp_mw is None:
            disp_mw = [[0.0 for _ in range(len(geom)*3)]]
            freq = [-1.0]
        else:
            assert len(freq) == len(disp_mw)

        if unit_xyz.lower() in ['bohr', 'au', 'a.u.']:
            pass
        elif unit_xyz.lower() in ['angstrom']:
            geom = [[g[0], (np.array(g[1]) * units.ANGSTROM).tolist()] for g in geom]
            disp_mw = [(np.array(d) * units.ANGSTROM).tolist() for d in disp_mw]
        else:
            raise NotImplementedError
        
        if unit_freq.lower() in ['cm1', 'cm', 'cm-1', 'kayser']:
            freq = (np.array(freq) * units.CM1).tolist()
        elif unit_freq.lower() in ['ev']:
            freq = (np.array(freq) * units.EV).tolist()
        elif unit_freq.lower() in ['au', 'a.u.', 'hartree']:
            pass
        else:
            raise NotImplementedError
        
        if unit_mass.lower() in ['amu', 'dalton']:
            disp_mw = [(np.array(d) * math.sqrt(units.DALTON)).tolist() for d in disp_mw]
        elif unit_mass.lower() in ['a.u.', 'au', 'emu']:
            pass
        else:
            raise NotImplementedError

        self.freq = freq
        self.disp_mw = disp_mw

        self.geom = geom
        self.atom = [a[0] for a in geom]
        self.coord = [c[1] for c in geom]
        self.bond = [[self.coord[i], self.coord[k]] for (i,k) in itertools.combinations(range(len(self.atom)), 2) \
                if (mendeleev.element(self.atom[i]).covalent_radius_pyykko +
                mendeleev.element(self.atom[k]).covalent_radius_pyykko)*1.0e-02 + 0.1> 
                dist(self.coord[i], self.coord[k]) / units.ANGSTROM]
        
        self.set_disp()

        self.nmode = len(self.disp_mw)
        self.natom = len(self.geom)

        self.set_disp()

        self.Q_mat = np.matrix(self.disp).T
    
    def set_disp(self):
        self.disp = []
        for d in self.disp_mw:
            self.disp.append(copy.deepcopy(d))
            for i, a in enumerate(self.atom):
                self.disp[-1][3*i] /= math.sqrt(mendeleev.element(a).atomic_weight * units.DALTON)
                self.disp[-1][3*i+1] /= math.sqrt(mendeleev.element(a).atomic_weight * units.DALTON)
                self.disp[-1][3*i+2] /= math.sqrt(mendeleev.element(a).atomic_weight * units.DALTON)

    def visualize(self, arrow_scale = 10, blender=False, atom_number=False):
        """Visualize vibrational modes.

        Vibration mode visualization with Matplotlib and tinker.

        Args:
            arrow_scale (float) : The scale of displacement arrows. Defaults to ``1``.
            blender (Optional[bool]) : View by blender.
            atom_number (Optional[bool]) : Plot atom number

        Examples:
            >>> sim = Vibration(geom, freq, disp)
            >>> sim.visualize()

        """
        if blender:
            raise NotImplementedError
        self.atom_number = atom_number

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
            command = lambda:inputData.plot(Canvas, ax, fig))
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

    def localize(self, option='Pipek-Mezy', window=400):
        r"""Localize vibrational modes.

        Vibration mode localization with local_cls module.

        - Pipek-Mezy metric
        .. math::

            \xi_{\mathrm{at}}\left(\widetilde{\boldsymbol{Q}}^{\mathrm{sub}}\right)=\sum_{p=1}^{k} \sum_{i=1}^{n}\left(\tilde{C}_{i p}^{\mathrm{sub}}\right)^{2}
            \\
            \tilde{C}_{i p}^{\mathrm{sub}}=\sum_{\alpha=x, y, z}\left(\tilde{Q}_{i \alpha, p}^{\mathrm{sub}}\right)^{2}

        - Boys metric
        .. math::

            \xi_{\mathrm{dist}}\left(\widetilde{\boldsymbol{Q}}^{\mathrm{sub}}\right)=\sum_{p=1}^{k}\left(\boldsymbol{R}_{p}^{\text {center }}\right)^{2}
            \\
            \boldsymbol{R}_{p}^{\text {center }}=\sum_{i=1}^{n} \tilde{C}_{i p}^{\mathrm{sub}} \boldsymbol{R}_{i}

        Args:
            option (str) : The metric of localization. ``Boys`` or ``Pipek-Mezy``.
            window (float) : window frequency in cm-1.

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

    def group_localize(self, mwhess : np.ndarray, 
        domain : List[List[float]],
        unit_xyz : Optional[str] ='au', 
        unit_omega : Optional[str] ='au',
        unit_mass : Optional[str] ='au'
        ) -> Tuple[List[List[float]], List[float]]:
        """Generate Group Localized Coordinate from mass-weighted hessian
        
        Args:
            mwhess (np.ndarray) : mass-weighted hessian
            domain (List[List[float]]): atomic domain such as [[0,1,2],[3,4]]
        """
        
        if unit_omega.lower() in ['cm1', 'cm', 'cm-1', 'kayser']:
            mwhess *= units.CM1**2
        elif unit_omega.lower() in ['ev']:
            mwhess *= units.EV**2
        elif unit_omega.lower() in ['au', 'a.u.', 'hartree']:
            pass
        else:
            raise NotImplementedError
        
        loc = GroupLocalizer(self, mwhess, domain)
        self.Q_mat, self.freq = loc.run()
        
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
        self.disp_mw = self.Q_mat.T.tolist()
        self.set_disp()
        return (self.disp_mw, self.freq)
