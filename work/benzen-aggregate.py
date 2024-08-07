#!/usr/bin/env python
# coding: utf-8

# # Quick Start

# System : H${}_2$CO
# 
# ## 1. Install and Satisfy requirements
# 
# See [README.md](https://github.com/KenHino/PyVibLocalizer/blob/main/README.md#installation).

# ## 2. Import modules

# In[1]:


import mendeleev
import numpy as np
import scipy


# In[2]:


try:
    import bpy
except ImportError:
    print('You cannot use Blender in Jupyter')
    print('You must select jupyter kernel as blender when you use blender')
    import matplotlib
try:
    import ase
    import ase.db
except ImportError:
    print('You cannot use ASE')


# In[3]:


try:
    import pyvib
except ImportError:
    print('Try execute in `src` directory or set $PYTHONPATH to pyvib')


# ## 3. Prepare Geometry
# 
# Define
# - atom element
# - coordinate
# 
# For example,
# ```python
# geom = [['C', (0.0, 0.0, 0.0)],
#         ['H', (1,0, 1.0, 1.0)]]
# ```
# 
# Here, I use precalculated `ase` database.

# In[4]:



geom = [
 ['C' ,(1.715018, 5.341687, 6.37253 )],
 ['H' ,(0.967674, 5.281308, 5.822208)],
 ['C' ,(2.713164, 4.351680, 6.318648)],
 ['H' ,(2.631569, 3.622774, 5.746772)],
 ['C' ,(7.237018, 5.341687, 6.37253 )],
 ['H' ,(6.489674, 5.281308, 5.822208)],
 ['C' ,(8.235164, 4.351680, 6.318648)],
 ['H' ,(8.153569, 3.622774, 5.746772)],
 ['C' ,(3.808113, 4.487670, 7.131918)],
 ['C' ,(9.330113, 4.487670, 7.131918)],
 ['H' ,(4.480340, 3.846341, 7.094559)],
 ['H' ,(10.00234, 3.846341, 7.094559)],
 ['C' ,(1.715018, 10.78129, 6.37253 )],
 ['H' ,(0.967674, 10.72091, 5.822208)],
 ['C' ,(2.713164, 9.791280, 6.318648)],
 ['H' ,(2.631569, 9.062374, 5.746772)],
 ['C' ,(7.237018, 10.78129, 6.37253 )],
 ['H' ,(6.489674, 10.72091, 5.822208)],
 ['C' ,(8.235164, 9.791280, 6.318648)],
 ['H' ,(8.153569, 9.062374, 5.746772)],
 ['C' ,(3.808113, 9.927270, 7.131918)],
 ['C' ,(9.330113, 9.927270, 7.131918)],
 ['H' ,(4.480340, 9.285941, 7.094559)],
 ['H' ,(10.00234, 9.285941, 7.094559)],
 ['C' ,(1.849344, 6.391530, 7.236809)],
 ['C' ,(1.849344, 11.83113, 7.236809)],
 ['H' ,(1.177116, 7.032859, 7.274168)],
 ['H' ,(1.177116, 12.47246, 7.274168)],
 ['C' ,(0.368383, 8.257313, 9.964712)],
 ['H' ,(-0.37896, 8.317692, 9.41439 )],
 ['C' ,(1.366528, 9.247320, 9.910829)],
 ['H' ,(1.284934, 9.976226, 9.338954)],
 ['C' ,(0.502708, 7.207470, 10.82899)],
 ['H' ,(-0.16952, 6.566141, 10.86635)],
 ['C' ,(7.371344, 6.391530, 7.236809)],
 ['C' ,(7.371344, 11.83113, 7.236809)],
 ['H' ,(6.699116, 7.032859, 7.274168)],
 ['H' ,(6.699116, 12.47246, 7.274168)],
 ['C' ,(2.595802, 8.061487, 11.58838)],
 ['C' ,(8.117802, 8.061487, 11.58838)],
 ['C' ,(3.942438, 5.537513, 7.996197)],
 ['C' ,(3.942438, 10.97711, 7.996197)],
 ['C' ,(9.464438, 5.537513, 7.996197)],
 ['C' ,(9.464438, 10.97711, 7.996197)],
 ['C' ,(5.890383, 8.257313, 9.964712)],
 ['H' ,(3.343147, 8.001108, 12.1387 )],
 ['H' ,(8.865147, 8.001108, 12.1387 )],
 ['H' ,(4.689783, 5.597892, 8.546519)],
 ['H' ,(4.689783, 11.03749, 8.546519)],
 ['H' ,(10.21178, 5.597892, 8.546519)],
 ['H' ,(10.21178, 11.03749, 8.546519)],
 ['H' ,(5.143038, 8.317692, 9.41439 )],
 ['C' ,(1.597657, 7.071480, 11.64226)],
 ['C' ,(7.119657, 7.071480, 11.64226)],
 ['C' ,(2.944292, 6.527520, 8.050079)],
 ['C' ,(2.944292, 11.96712, 8.050079)],
 ['C' ,(8.466292, 6.527520, 8.050079)],
 ['C' ,(8.466292, 11.96712, 8.050079)],
 ['C' ,(6.888528, 9.247320, 9.910829)],
 ['H' ,(1.679251, 6.342574, 12.21414)],
 ['H' ,(7.201251, 6.342574, 12.21414)],
 ['H' ,(3.025887, 7.256426, 8.621955)],
 ['H' ,(3.025887, 12.69603, 8.621955)],
 ['H' ,(8.547887, 7.256426, 8.621955)],
 ['H' ,(8.547887, 12.69603, 8.621955)],
 ['H' ,(6.806934, 9.976226, 9.338954)],
 ['C' ,(2.461477, 9.111330, 10.7241 )],
 ['C' ,(7.983477, 9.111330, 10.7241 )],
 ['C' ,(6.024708, 7.207470, 10.82899)],
 ['H' ,(3.133704, 9.752659, 10.68674)],
 ['H' ,(8.655704, 9.752659, 10.68674)],
 ['H' ,(5.352481, 6.566141, 10.86635)]]


# In[5]:




# ## 4. Prepare mass-weighted hessian
# 
# Define
# - mass-weighted hessian (unit is a.u.)
# 
# Mass-weighted hessian is matrix $M_{ij}=\frac{\partial^2 E}{\partial\sqrt{m_i}x_i\partial\sqrt{m_j}x_j}$ of the second derivative of energy $E$ in terms of mass-weighted coordinates $\sqrt{m_i}x_i$. You can also obtain mass-weighted hessian from harmonic frequency and displacement vectors.
# 
# You can use `pyvib.read_fchk_g16()` for [Gaussian16](https://gaussian.com/gaussian16/) or `pyvib.read_minfo()` for [SINDO](https://tms.riken.jp/research/software/sindo/)

# In[6]:




# In[7]:




# ## 5. Set PyViblocalizer
# - Set geometry and (hessian or (displacement vector and frequency))
# - Input units can be specified in options, such as `unit_mass='AMU'` 
# 
# In large system, it may takes a few minutes

# In[8]:


vib = pyvib.Vibration(geom, unit_xyz='angstrom')


# ## 6.1. Visualize in Tkinter+matplotlib

# ## 6.2. Visualize in Blender

# In[9]:


try:
    import bpy
    vib.visualize(blender=True)
except ImportError:
    '''Not recommended'''
    vib.visualize(blender=False)


