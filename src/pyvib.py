import itertools
from collections import defaultdict, Counter
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
                line = line.replace(',','').replace('E','e')
                words = line.split()
                if len(words) == 6:
                    geom.append([words[0], [float(w) for w in words[3:6]]])
                else:
                    natom = int(words[0])
                continue

        elif read_freq:
            if line[-2] == ' ':
                line = line.replace(',','').replace('e','e')
                words = line.split()
                freq.extend([float(w) for w in words])
                continue

        elif read_disp:
            if line[-2] == ' ':
                line = line.replace(',','').replace('e','e')
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
                if (data[self.atom[i]][1]+data[self.atom[k]][1]) + 0.15 > dist(self.coord[i], self.coord[k])]

        
        self.nmode = len(self.disp)
        self.natom = len(self.geom)

        self.Q_mat = np.matrix(self.disp).T


    def visualize(self, arrow_scale = 10):
        ### root object ###
        root = tk.Tk()
        root.title("PyVibVisualizer")

        inputFrame = ttk.Frame(root)
        buttonFrame = ttk.Frame(root)
        graphFrame = ttk.Frame(root)

        ### inputFrame ###
        inputData = visualizer(inputFrame,self, arrow_scale)
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

    def localize(self, *, option='Pipek-Mezy', window=400):
        loc = localizer(self, option, window)
        self.Q_mat, self.freq = loc.run()
        for x in range(self.nmode):
            for y in range(self.natom*3):
                self.disp[x][y] = self.Q_mat[y,x]

        print(option,'metric = ',loc.metric(self.Q_mat))

