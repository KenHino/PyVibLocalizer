# PyVibLocalizer
## create and visualize local mode from normal mode!

J. Chem. Phys. 130, 084106 (2009); https://doi.org/10.1063/1.3077690

To maximize, metrices, this program used `scipy.optimized.minimize_scalar`, in jacobi sweep method. The choice of unitary rotation mode pair is based on window frequency.: J. Chem. Phys. 145, 124112 (2016); https://doi.org/10.1063/1.4963109 Default window frequency set as 400 cm-1.


- coordinate unit must be Bohr

- `.minfo` file is made by SINDO(https://tms.riken.jp/research/software/sindo/)

- to run this program, execute `python {$mol_name}.py` (example; moL_name=ch2o, c6h8)
