import itertools
import math

import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from atom_data import data
from visual_cls import visualizer
from local_cls import localizer

def dist(A, B):
    dum = 0.0
    for a,b in zip(A,B):
        dum += (a-b)**2
    return math.sqrt(dum)

def read_minfo(file_name):
    geom = []
    freq = []
    disp = []
    return geom, freq, disp

def write_minfo(file_name, geom, freq, disp):
    a =0

def read_gout(file_name):
    geom = []
    freq = []
    disp = []
    return geom, freq, disp

def write_gout(file_name, geom, freq, disp):
    a=0



class simulator:
    def __init__(self, geom, freq, disp, unit='Bohr'):
        self.unit = unit
        if unit != 'Bohr':
            assert False, 'Not yet implemented'
        self.freq = freq
        self.disp = disp
        self.geom = geom
        self.atom = [a[0] for a in geom]
        self.coord = [c[1] for c in geom]
        self.bond = [[self.coord[i], self.coord[k]] for (i,k) in itertools.combinations(range(len(self.atom)), 2) \
                if (data[self.atom[i]][1]+data[self.atom[k]][1]) > dist(self.coord[i], self.coord[k])]

        
        self.nmode = len(self.disp)
        self.natom = len(self.geom)

        self.Q_mat = np.matrix(self.disp).T


    def visualize(self):
        ### root object ###
        root = tk.Tk()
        root.title("PyVibVisualizer")

        inputFrame = ttk.Frame(root)
        buttonFrame = ttk.Frame(root)
        graphFrame = ttk.Frame(root)

        ### inputFrame ###
        inputData = visualizer(inputFrame,self)
        inputFrame.pack()

        ### buttonFrame ###
        ButtonWidth = 10
        UpdateButton = tk.Button(buttonFrame, text="View", width=ButtonWidth, \
            command = lambda:inputData.plot(Canvas, ax, fig))
        UpdateButton.grid(row = 0, column = 0)
        buttonFrame.pack()

        #graph initialize
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect((1,1,1))
        ax.axis('off')
        Canvas = FigureCanvasTkAgg(fig, master = graphFrame)
        Canvas.get_tk_widget().pack()
        graphFrame.pack()

        #continue
        root.mainloop()

    def localize(self, *, option='Pipek-Mezy'):
        loc = localizer(self, option)
        self.Q_mat, self.freq = loc.run()
        for x in range(self.nmode):
            for y in range(self.natom*3):
                self.disp[x][y] = self.Q_mat[y,x]

        print(option,'metric = ',loc.metric(self.Q_mat))

