from scipy.constants import physical_constants as pc
"""  Units 
Example
   1 [angstrom] = ANGSTROM     [a.u.]
   1 [a.u.]     = 1 / ANGSTROM [angstrom]
"""

ANGSTROM = 1 / pc['Bohr radius'][0] * 1.0e-10
DALTON = pc['atomic mass constant'][0] / pc['electron mass'][0]
EV =  1 / pc['Hartree energy in eV'][0] 
FS = 1 / pc['atomic unit of time'][0] * 1.0e-15 
CM1 = 1 / (pc['atomic unit of energy'][0] /(pc['speed of light in vacuum'][0] * 1.e02) / pc['Planck constant'][0] )

if __name__ == '__main__':
    print(f'ANGSTROM : {ANGSTROM} [angstrom / a.u.]')
    print(f'DALTON   : {DALTON} [AMU / a.u.]')
    print(f'EV       : {EV} [eV / a.u.]')
    print(f'FS       : {FS} [fs / a.u.]')
    print(f'CM1      : {CM1} [cm^-1 / a.u.]')
