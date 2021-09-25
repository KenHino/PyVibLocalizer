# PyVibLocalizer
## create and visualize local mode from normal mode!

J. Chem. Phys. 130, 084106 (2009); https://doi.org/10.1063/1.3077690

- Pipek-Mezy metric
    \
    <img src="https://latex.codecogs.com/gif.latex?\space\xi_{\mathrm{at}}\left(\widetilde{\boldsymbol{Q}}^{\mathrm{sub}}\right)=\sum_{p=1}^{k} \sum_{i=1}^{n}\left(\tilde{C}_{i p}^{\mathrm{sub}}\right)^{2}" />
    \
    <img src="https://latex.codecogs.com/gif.latex?\space\tilde{C}_{i p}^{\mathrm{sub}}=\sum_{\alpha=x, y, z}\left(\tilde{Q}_{i \alpha, p}^{\mathrm{sub}}\right)^{2}" />
- Boys metric
    \
    <img src="https://latex.codecogs.com/gif.latex?\space\xi_{\mathrm{dist}}\left(\widetilde{\boldsymbol{Q}}^{\mathrm{sub}}\right)=\sum_{p=1}^{k}\left(\boldsymbol{R}_{p}^{\text {center }}\right)^{2}" />
    \
    <img src="https://latex.codecogs.com/gif.latex?\space\boldsymbol{R}_{p}^{\text {center }}=\sum_{i=1}^{n} \tilde{C}_{i p}^{\mathrm{sub}} \boldsymbol{R}_{i}"/>

to maximize, metrices, this program used `scipy.optimized.miniimize`, not jacobi sweep method.


- coordinate unit must be Bohr

- `.minfo` file is made by SINDO(https://tms.riken.jp/research/software/sindo/)

- to run this program, execute `python main.py`