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


geom, _ = pyvib.read_fchk_g16('./sample/ch2o.fchk')


# In[5]:


geom


# ## 4. Prepare mass-weighted hessian
# 
# Define
# - mass-weighted hessian (unit is a.u.)
# 
# Mass-weighted hessian is matrix $M_{ij}=\frac{\partial^2 E}{\partial\sqrt{m_i}x_i\partial\sqrt{m_j}x_j}$ of the second derivative of energy $E$ in terms of mass-weighted coordinates $\sqrt{m_i}x_i$. You can also obtain mass-weighted hessian from harmonic frequency and displacement vectors.
# 
# You can use `pyvib.read_fchk_g16()` for [Gaussian16](https://gaussian.com/gaussian16/) or `pyvib.read_minfo()` for [SINDO](https://tms.riken.jp/research/software/sindo/)

# In[6]:


_, mw_hess = pyvib.read_fchk_g16('./sample/ch2o.fchk')


# In[7]:


mw_hess.shape


# ## 5. Set PyViblocalizer
# - Set geometry and (hessian or (displacement vector and frequency))
# - Input units can be specified in options, such as `unit_mass='AMU'` 
# 
# In large system, it may takes a few minutes

# In[8]:


vib = pyvib.Vibration(geom, mw_hess=np.array(mw_hess), 
                      unit_omega='hartree', unit_mass='a.u.')


# ## 6.1. Visualize in Tkinter+matplotlib

# ## 6.2. Visualize in Blender

# In[9]:


try:
    import bpy
    vib.visualize(blender=True)
except ImportError:
    '''Not recommended'''
    vib.visualize(blender=False)


# ![image.png](attachment:8ca6e183-f7ec-4df6-b843-1a76502e4f8d.png)

# ## 7. Localization
# Returns displacement vecotors in a.u. and frequency in a.u.

# ## 7.1 Group Localization

# In[10]:


vib = pyvib.Vibration(geom, mw_hess=np.array(mw_hess), 
                      unit_omega='hartree', unit_mass='a.u.')
disp, freq = vib.group_localize(domain=[[0,1],[2,3]], mw_hess=np.array(mw_hess), 
                                   unit_omega='hartree', unit_mass='a.u.')


# In[11]:


try:
    import bpy
    vib.visualize(blender=True)
except ImportError:
    vib.visualize(blender=False)


# ![image.png](attachment:8f0c42d9-04f9-4bd7-b1ca-4f6953028863.png)

# ## 7.2. Metric Localization

# In[12]:


vib = pyvib.Vibration(geom, mw_hess=np.array(mw_hess), 
                      unit_omega='hartree', unit_mass='a.u.')
disp, freq = vib.localize(option='Boys', window= 500)


# In[13]:


try:
    import bpy
    vib.visualize(blender=True)
except ImportError:
    vib.visualize(blender=False)


# ![image.png](attachment:baa00b52-c7eb-40c8-8d85-69da6f45250b.png)
