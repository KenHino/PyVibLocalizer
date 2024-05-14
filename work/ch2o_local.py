import numpy as np
from pyvib import Vibration, read_minfo #geometry, freq, displacement vector

geom = [['O',     [   1.13890361,      -0.00029065,      -0.00005828]],
        ['C',     [  -1.13628333,       0.00028780,       0.00005557]],
        ['H',     [  -2.27284653,       0.46535179,      -1.72435469]],
        ['H',     [  -2.27284922,      -0.46416573,       1.72461804]]]#'Bohr'

'''NOT dimensionless coordinate frequency'''
freq = [ 1.18646597e+03,  1.25305305e+03,  1.51535414e+03,  1.83295234e+03,  2.86496842e+03, 2.91731717e+03] # cm-1

'''Normalized mass-weighted displacement vector'''
disp = [[-3.77832663e-05, -1.44206524e-01, -3.90016272e-02,
          1.37330386e-04,  4.99831945e-01,  1.34928755e-01,
          1.07701428e-04, -5.75125668e-01, -1.55109452e-01,
         -4.31056129e-04, -5.75116981e-01, -1.55104398e-01],
        [ 4.02853309e-06, -7.18872634e-02,  2.66462060e-01,
          1.15371251e-06,  1.15979195e-01, -4.29385874e-01,
         -5.61828755e-01, -5.67706738e-02,  2.10097471e-01,
          5.61808725e-01, -5.70452510e-02,  2.10019470e-01],
        [-3.00028578e-01,  7.61788379e-05,  1.97139367e-05,
          8.67375363e-03,  1.06048618e-06, -3.88969539e-06,
          5.82655326e-01,  8.82669590e-02, -3.28116731e-01,
          5.82671009e-01, -8.85741002e-02,  3.28051616e-01],
        [ 6.13995166e-01, -1.55619827e-04, -3.30190677e-05,
         -7.48945014e-01,  1.90049455e-04,  4.33612467e-05,
          6.91488917e-02,  4.21543319e-02, -1.56489074e-01,
          6.91423931e-02, -4.21901622e-02,  1.56470992e-01],
        [-4.07764541e-03,  7.56488909e-07,  2.30879745e-06,
         -1.97855754e-01,  8.38679699e-05, -1.04146966e-04,
          3.49610346e-01, -1.55904927e-01,  5.78117351e-01,
          3.49360671e-01,  1.55612517e-01, -5.77767176e-01],
        [ 2.98490039e-06, -8.90723115e-04,  3.30468722e-03,
          7.19831690e-05,  8.64782849e-02, -3.20876399e-01,
          3.51718386e-01, -1.47466018e-01,  5.46796360e-01,
         -3.51978664e-01, -1.47389841e-01,  5.47263255e-01]]#'Bohr(mass(AMU???) weighted)'

U = np.load('U.npy')
disp = U @ np.array(disp)
disp[3:5,:] *= -1

freq = np.diag(U@np.diag(freq)@U.T)

#geom, freq, disp = read_minfo("sample/ch2o.minfo")


#sim = Vibration(geom, freq, disp)

#sim.localize(option='Boys')
#sim.localize(option='Pipek-Mezy', window=500)
#sim.visualize()
import numpy as np
#mw_hess = np.array(disp).T@np.diag(np.array(freq)**2)@np.array(disp)
#sim = Vibration(geom, mw_hess=np.array(mw_hess), unit_omega='cm-1', unit_mass='AMU')
sim = Vibration(geom, freq=freq, mw_disp=disp, unit_omega='cm-1')
#sim.group_localize(domain=[[0,1,2,3]], mw_hess=mw_hess, unit_omega='cm-1', unit_mass='AMU')
sim.visualize(blender=True)
