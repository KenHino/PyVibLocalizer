# PyVibLocalizer
## create and visualize local mode from normal mode!

J. Chem. Phys. 130, 084106 (2009); https://doi.org/10.1063/1.3077690

- Pipek-Mezy metric
    $$\xi_{\mathrm{at}}\left(\widetilde{\boldsymbol{Q}}^{\mathrm{sub}}\right)=\sum_{p=1}^{k} \sum_{i=1}^{n}\left(\tilde{C}_{i p}^{\mathrm{sub}}\right)^{2}$$
    $$\tilde{C}_{i p}^{\mathrm{sub}}=\sum_{\alpha=x, y, z}\left(\tilde{Q}_{i \alpha, p}^{\mathrm{sub}}\right)^{2}$$
- Boys metric
    $$\xi_{\mathrm{dist}}\left(\widetilde{\boldsymbol{Q}}^{\mathrm{sub}}\right)=\sum_{p=1}^{k}\left(\boldsymbol{R}_{p}^{\text {center }}\right)^{2}$$
    $$\boldsymbol{R}_{p}^{\text {center }}=\sum_{i=1}^{n} \tilde{C}_{i p}^{\mathrm{sub}} \boldsymbol{R}_{i}$$

To maximize, metrices, this program used `scipy.optimized.minimize_scalar`, in jacobi sweep method. The choice of unitary rotation mode pair is based on window frequency.: J. Chem. Phys. 145, 124112 (2016); https://doi.org/10.1063/1.4963109 Default window frequency set as 400 cm-1.


- coordinate unit must be Bohr

- `.minfo` file is made by SINDO(https://tms.riken.jp/research/software/sindo/)

- to run this program, execute `python {$mol_name}.py` (example; moL_name=ch2o, c6h8)