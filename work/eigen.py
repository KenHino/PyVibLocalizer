#!/usr/bin/env python
# coding: utf-8

# # Eigen Cation

# ## 1. Install and Satisfy requirements
# 
# See [documenstaion](https://kenhino.github.io/PyVibLocalizer/README.html).

# ## 2. Import modules

# In[14]:


import mendeleev
import numpy as np
import scipy


# In[15]:


try:
    import bpy
except ImportError:
    print('You cannot use Blender in Jupyter')
    print('You must select jupyter kernel as blender when you use blender')
    import mendeleev
import matplotlib
try:
    import ase
    import ase.db
except ImportError:
    print('You cannot use ASE')


# In[16]:


try:
    import pyvib
except ImportError:
    print('Try execute in `src` directory or set $PYTHONPATH to pyvib')


# ## 3. Prepare Geometry
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
# Here, I use pre-calculated coordinate.

# In[17]:


geom, mw_hess = pyvib.read_fchk_g16('./sample/eigen.fchk')


# In[18]:


geom


# ## 4. Prepare mass-weighted hessian
# 
# Define
# - mass-weighted hessian (unit is a.u.)
# 
# Mass-weighted hessian is matrix $M_{ij}=\frac{\partial^2 E}{\partial\sqrt{m_i}x_i\partial\sqrt{m_j}x_j}$ of the second derivative of energy $E$ in terms of mw coordinates $\sqrt{m_i}x_i$. You can also obtain mw-hessian from harmonic frequency and displacement vectors.
# 
# You can use `pyvib.read_minfo()` for [SINDO](https://tms.riken.jp/research/software/sindo/) or `pyvib.read_fchk_g16()` for [Gaussian16](https://gaussian.com/gaussian16/)

# In[19]:


mw_hess


# In[20]:


mw_hess.shape


# ## 5. Set PyViblocalizer
# - Set geometry and (hessian or (displacement vector and frequency))
# - Input units can be specified in options, such as `unit_mass='AMU'` 
# 
# In large system, it may takes a few minutes

# In[21]:


vib = pyvib.Vibration(geom, mw_hess=np.array(mw_hess), 
                      unit_xyz='bohr', unit_mass='a.u.')
import units
np.array(vib.freq) / units.CM1


# ## 7. Localization

# ## 7.1 Group Localization

# In[22]:


#disp, freq = vib.group_localize(domain=[[0,1,2],[3,4,5],[6,7,8,9],[10,11,12]], mw_hess=np.array(mw_hess), unit_mass='a.u.')
disp, freq = vib.group_localize(domain=[[0,1,2],[3,4,5],[7,8,9],[10,11,12]], mw_hess=np.array(mw_hess), unit_mass='a.u.')
#disp, freq = vib.group_localize(domain=[[0,1,2,3,4,5,6,7,8,9,10,11,12]], mw_hess=np.array(mw_hess), unit_mass='a.u.')
print(np.array(disp))
print(np.array(freq) / units.CM1)
np.save('freq_au_local', np.array(freq))
np.save('disp_au_local', np.array(disp))


# In[23]:


try:
    import bpy
    vib.visualize(blender=True, arrow_scale=150.0)
except ImportError:
    vib.visualize(blender=False)


# ![image.png](attachment:562076bc-5e32-4c13-99b6-ba47805d57d7.png)
