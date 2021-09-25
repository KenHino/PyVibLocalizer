import numpy as np
from scipy.optimize import minimize
import math

class localizer:
    def __init__(self, simu, option):
        self.geom = simu.geom
        self.Q_mat = np.matrix(simu.disp).T
        self.R_mat =  np.array(simu.coord)
        self.nmode = len(simu.disp)
        self.natom = len(self.geom)
        self.option = option
        
        self.freq = simu.freq
        self.hamiltonian = np.diag(-np.power(np.array(self.freq),2))


    def metric(self, Q_mat):
        if self.option == 'Pipek-Mezy':
            dum = 0.0
            for pmode in range(self.nmode):
                for iatom in range(self.natom):
                    dum += np.sum(np.power(Q_mat[3*iatom:3*(iatom+1), pmode],2))**2

        elif self.option == 'Boys':
            dum = 0.0
            for pmode in range(self.nmode):
                R_center = np.zeros(3)
                for iatom in range(self.natom):
                    R_center += np.sum(np.power(Q_mat[3*iatom:3*(iatom+1), pmode],2))*self.R_mat[iatom]
                dum += np.linalg.norm(R_center, ord=2)
        else:
            assert False

        return dum

    def construct_unitary(self, pair_theta_dic):
        identity = np.eye(self.nmode)
        for pair, theta in pair_theta_dic.items():
            pmode = pair[0]
            qmode = pair[1]
            rotate_mat = np.eye(self.nmode)
            rotate_mat[pmode][pmode] = rotate_mat[qmode][qmode] = np.cos(theta)
            rotate_mat[pmode][qmode] = np.sin(theta)
            rotate_mat[qmode][pmode] = - np.sin(theta)
            identity = np.dot(identity, rotate_mat)
        return identity


    def run(self):

        def objective(array):
            pair_theta = {}
            for pmode in range(self.nmode-1):
                pair_theta[(pmode, pmode+1)] = array[pmode]
            rotate_matrix = self.construct_unitary(pair_theta)
            Q_mat = np.dot(self.Q_mat, rotate_matrix)
            return - self.metric(Q_mat)

        result = minimize(objective, np.zeros(self.nmode-1))
        print(result.message)
        if not result.success:
            assert False

        pair_theta = {}
        for pmode in range(self.nmode-1):
            pair_theta[(pmode, pmode+1)] = result.x[pmode]

        rotate_matrix = self.construct_unitary(pair_theta)

        self.Q_mat = np.dot(self.Q_mat, rotate_matrix)

        np.set_printoptions(formatter={'float': '{:>8.0f}'.format})
        print('\n','localized hessian [cm-2]')
        hamiltonian = np.dot(np.dot(rotate_matrix,  self.hamiltonian), rotate_matrix.T)
        print(hamiltonian,'\n')

        return self.Q_mat, list(np.sqrt(-np.diag(hamiltonian)))


