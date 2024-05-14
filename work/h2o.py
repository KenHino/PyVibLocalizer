from pyvib import Vibration, read_minfo #geometry, freq, displacement vector

geom = [['O',     [    0.000000,      0.000000,      0.236470]],
        ['H',     [    0.000018,      1.429721,     -0.906721]],
        ['H',     [   -0.000018,     -1.429721,     -0.906721]]]#'Bohr'

'''NOT dimensionless coordinate frequency'''
freq = [ 1658.83, 3750.33, 3851.95] # cm-1

'''cartesian displacement vector'''

disp = [[0.00000,  0.00000, -0.06761, 0.00001, 0.41458, 0.53657, -0.00001, -0.41458, 0.53657],
        [0.00000,  0.00000, 0.04923, 0.00001, 0.56937, -0.39069, -0.00001, -0.56937, -0.39069],
        [0.00000, -0.06679, 0.00000, 0.00001, 0.53009, -0.42386, 0.00001, 0.53009, 0.42386]]

#geom, freq, disp = read_minfo("sample/ch2o.minfo")


sim = Vibration(geom, freq, disp=disp, unit_xyz='bohr', unit_omega='cm-1')

#sim.localize(option='Boys')
#sim.localize(option='Pipek-Mezy', window=500)
#sim.visualize()
sim.visualize(blender=True, arrow_scale=2)
