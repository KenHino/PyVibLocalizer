import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import math
import mendeleev

from mpl_toolkits.mplot3d import Axes3D

import units

class Visualizer(ttk.Frame):
    def __init__(self, inputFrame, vibration, arrow_scale):

        self.freq = vibration.freq
        self.disp = vibration.disp
        self.geom = vibration.geom
        self.bond = vibration.bond
        self.coord = vibration.coord
        self.atom = vibration.atom
        self.atom_number = vibration.atom_number
        self.natom = len(self.geom)
        self.scale = arrow_scale

        self.trgFunc = tk.StringVar()
        self.ComboboxTrgFunc = ttk.Combobox(inputFrame, 
                        textvariable=self.trgFunc, width=13)
        select = (f"{k+1} : "+\
            f"{int(self.freq[k] / units.CM1) if not math.isnan(self.freq[k]) else 0.0}cm-1" 
                    for k in range(len(self.freq)))
        self.ComboboxTrgFunc['values'] = tuple(select)
        self.ComboboxTrgFunc.insert(tk.END, "select mode")
        self.ComboboxTrgFunc.pack(side=tk.LEFT)

    def plot_sphere(self, r, center, color, ax):
        # Make data
        _N = 30
        u = np.linspace(0, 2 * np.pi, _N)
        v = np.linspace(0, np.pi, _N)
        x = r * np.outer(np.cos(u), np.sin(v)) + center[0]
        y = r * np.outer(np.sin(u), np.sin(v)) + center[1]
        z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]

        # Plot the surface
        ax.plot_surface(x, y, z, color=color, rcount=_N, 
                        ccount=_N, antialiased=False, alpha=0.3)

    def plot_bond(self, bond, ax):
        x = np.linspace(bond[0][0] / units.ANGSTROM, bond[1][0] / units.ANGSTROM, 30)
        y = np.linspace(bond[0][1] / units.ANGSTROM, bond[1][1] / units.ANGSTROM, 30)
        z = np.linspace(bond[0][2] / units.ANGSTROM, bond[1][2] / units.ANGSTROM, 30)
        ax.plot(x, y, z, ms=2, linewidth=4, color='gray')

    def plot_arrow(self, starts, vectors, ax):
        start = np.array(starts) / units.ANGSTROM
        vec = np.array(vectors).reshape(self.natom, 3) / units.ANGSTROM
        for xyz, uvw in zip(start, vec):
            norm = np.linalg.norm(uvw) * self.scale
            ax.quiver(*xyz,*uvw, pivot = 'tail', 
                length = norm, linewidths=4, color='green')
    
    def plot_number(self, number, coord, ax):
        ax.text(*coord, number, size=15)

    def plot(self, canvas, ax, fig):
        ax.cla()
        fig.clf()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect((1,1,1))
        ax.axis('off')
        for i, atom in enumerate(self.geom):
            self.plot_sphere(mendeleev.element(atom[0]).atomic_radius * 1.0e-02, 
                (np.array(atom[1]) / units.ANGSTROM).tolist(), 
                mendeleev.element(atom[0]).molcas_gv_color, ax)
            if self.atom_number:
                self.plot_number(i, (np.array(atom[1]) / units.ANGSTROM).tolist(), ax)
        for bond in self.bond:
            self.plot_bond(bond, ax)

        index = int(self.ComboboxTrgFunc.get().split()[0]) - 1
        self.plot_arrow(self.coord, self.disp[index], ax)
        
        margin = 0.5
        maximum = -1.e10
        minimum =  1.e10
        for atom in self.geom:
            maximum = max(maximum, max(atom[1]))
            minimum = min(minimum, min(atom[1]))

        ax.set_xlim(minimum - margin, maximum + margin)
        ax.set_ylim(minimum - margin, maximum + margin)
        ax.set_zlim(minimum - margin, maximum + margin)
        ax.set_xlabel('$\AA$')
        ax.set_ylabel('$\AA$')
        ax.set_zlabel('$\AA$')
        
        canvas.draw()
