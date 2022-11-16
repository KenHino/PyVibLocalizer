# PyVibLocalizer
**visualize atomic mode and genetrate displacement vectors (normal/local mode, arbitary path)**

<img src="./_docs/pic/ch2o.png" width="400">

- [Documentation](https://kenhino.github.io/PyVibLocalizer/index.html#)

## Installation
```bash
$ git clone https://github.com/KenHino/PyVibLocalizer
```

- Requirements
    - must
        ```bash
        $ pip install mendeleev numpy scipy
        ```

    - better
        - [Blender](https://www.blender.org/)

            If WSL2 or Ubuntu,
            ```bash
            $ sudo apt install blender # If Mac OS brew install --cask blender
            ```
            In your `.bashrc` files
            ```
            set PYTHONHOME=`which python3` # May be not not required
            alias blender='/usr/bin/blender --python-use-system-env'
            ```
            where `/usr/bin/blender` is your blender installed PATH. In Mac OS, this may be `/Applications/Blender.app/Contents/MacOS/Blender`.
            
            When you execute `main.py`
            ```
            $ blender --python main.py
            ```
            When install package in blender-python in Mac OS
            ```
            $ <BPYTHON> -m pip install <PACKAGE>
            ```
            where \<BPYTHON\> may be `/Applications/Blender.app/Contents/Resources/3.3/python/bin/python3.10`, 
            \<PACKAGE\>=`numpy, scipy, mendeleev, ase`.
            
            When you use Jupyter Notebook
            ```bash
            $ pip install blender-notebook
            $ blender_notebook install --blender-exec="/usr/bin/blender"
            ```
            where `/usr/bin/blender` is your blender installed PATH. In Mac OS, this may be `/Applications/Blender.app/Contents/MacOS/Blender`.
            
            If required,
            ```bash
            $ sudo apt install subversion # If Mac OS brew install svn
            $ pip install future_fstrings
            $ pip install bpy
            ```
        - Matplotlib+tkinter
            ```
            $ pip install matplotlib
            ```
            Tf you use WSL2, below command may be needed.
            ```bash
            $ sudo apt-get install python3-tk
            ```


        - [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
            ```bash
            $ pip install ase
            ```


## Local mode
[J. Chem. Phys. 130, 084106 (2009)](https://doi.org/10.1063/1.3077690)

- Pipek-Mezy metric

$$\xi_{\mathrm{at}}\left(\widetilde{\boldsymbol{Q}}^{\mathrm{sub}}\right)=\sum_{p=1}^{k} \sum_{i=1}^{n}\left(\tilde{C}_{i p}^{\mathrm{sub}}\right)^{2}$$

$$\tilde{C}_{i p}^{\mathrm{sub}}=\sum_{\alpha=x, y, z}\left(\tilde{Q}_{i \alpha, p}^{\mathrm{sub}}\right)^{2}$$

- Boys metric

$$\xi_{\mathrm{dist}}\left(\widetilde{\boldsymbol{Q}}^{\mathrm{sub}}\right)=\sum_{p=1}^{k}\left(\boldsymbol{R}_{p}^{\text {center }}\right)^{2}$$

$$\boldsymbol{R}_{p}^{\text {center }}=\sum_{i=1}^{n} \tilde{C}_{i p}^{\mathrm{sub}} \boldsymbol{R}_{i}$$

To maximize, metrices, this program used `scipy.optimized.minimize_scalar`, in jacobi sweep method. The choice of unitary rotation mode pair is based on window frequency.: [J. Chem. Phys. 145, 124112 (2016)](https://doi.org/10.1063/1.4963109) Default window frequency set as 400 cm-1.

## Group Localized Coordinate
Group localized coordinate is one of the local mode, which diagonalize subspace mass-weighted hessian whose subspace is defined by your selction of atoms.


## Quick Start
- See [here](https://kenhino.github.io/PyVibLocalizer/quick_start.html)
- To run this program, execute `python {$mol_name}.py` or `blender --python {$mol_name}.py` (example; mol_name=src/ch2o) 
