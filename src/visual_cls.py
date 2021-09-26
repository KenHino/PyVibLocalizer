import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

from mpl_toolkits.mplot3d import Axes3D

from atom_data import data

class visualizer(ttk.Frame):
    def __init__(self, inputFrame, simulator):

        self.freq = simulator.freq
        self.disp = simulator.disp
        self.geom = simulator.geom
        self.bond = simulator.bond
        self.coord = simulator.coord
        self.atom = simulator.atom
        self.natom = len(self.geom)


        self.trgFunc = tk.StringVar()
        self.ComboboxTrgFunc = ttk.Combobox(inputFrame, textvariable=self.trgFunc, width=13)
        select = ("{} : {}cm-1".format(k+1, int(self.freq[k])) for k in range(len(self.freq)))
        self.ComboboxTrgFunc['values'] = tuple(select)
        self.ComboboxTrgFunc.insert(tk.END, "select mode")
        self.ComboboxTrgFunc.pack(side=tk.LEFT)

    def plot_sphere(self,r, center, color, ax):
        # Make data
        _N = 30
        u = np.linspace(0, 2 * np.pi, _N)
        v = np.linspace(0, np.pi, _N)
        x = r * np.outer(np.cos(u), np.sin(v)) + center[0]
        y = r * np.outer(np.sin(u), np.sin(v)) + center[1]
        z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]

        # Plot the surface
        ax.plot_surface(x, y, z,color=color,rcount=_N, ccount=_N, antialiased=False)

    def plot_bond(self, bond, ax):
        x = np.linspace(bond[0][0], bond[1][0], 30)
        y = np.linspace(bond[0][1], bond[1][1], 30)
        z = np.linspace(bond[0][2], bond[1][2], 30)
        ax.plot(x,y,z,ms=2,linewidth=4, color='gray')

    def plot_arrow(self, starts, vector, ax, scale= 10):
        star = np.array(starts).T
        vec = np.array(vector).reshape(self.natom,3).T
        for x,y,z,u,v,w in zip(*star, *vec):
            norm = pow(u*u+v*v+w*w, 0.5)*scale
            ax.quiver(x,y,z,u,v,w, pivot = 'tail', length = norm, linewidths=4, color='green')


    def plot(self, canvas, ax, fig):
        ax.cla()
        fig.clf()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect((1,1,1))
        ax.axis('off')
        for atom in self.geom:
            self.plot_sphere(0.9,atom[1], data[atom[0]][2], ax)

        for bond in self.bond:
            self.plot_bond(bond,ax)

        index = int(self.ComboboxTrgFunc.get().split()[0]) - 1
        self.plot_arrow([coord[1] for coord in self.geom], self.disp[index], ax)



        margin = 3
        maximum = -1.e10
        minimum =  1.e10
        for atom in self.geom:
            maximum = max(maximum, max(atom[1]))
            minimum = min(minimum, min(atom[1]))


        ax.set_xlim(minimum-margin,maximum+margin)
        ax.set_ylim(minimum-margin,maximum+margin)
        ax.set_zlim(minimum-margin,maximum+margin)
        ax.set_xlabel('Bohr')
        ax.set_ylabel('Bohr')
        ax.set_zlabel('Bohr')
        
        canvas.draw()
