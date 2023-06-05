import itertools
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.linalg import eigh
from typing import List, Optional

class Localizer:
    def __init__(self, vib, option: str, window: float):
        self.geom = vib.geom
        self.Q_mat = np.matrix(vib.disp).T
        self.R_mat = np.array(vib.coord)
        self.nmode = len(vib.disp)
        self.natom = len(self.geom)
        self.option = option
        
        self.freq = vib.freq
        self.window = window
        self.unitary = np.eye(self.nmode)
        self.hess = np.diag(np.power(np.array(self.freq), 2))

    def metric(self, Q_mat: np.ndarray, check_mode: Optional[List[float]] ='all'):
        if check_mode == 'all':
            check_mode = [pmode for pmode in range(self.nmode)]
        if self.option.lower() == 'pipek-mezey':
            dum = 0.0
            for pmode in check_mode:
                for iatom in range(self.natom):
                    dum += np.sum(np.power(
                        Q_mat[3*iatom:3*(iatom+1), pmode], 2))**2

        elif self.option.lower() == 'boys':
            dum = 0.0
            for pmode in check_mode:
                R_center = np.zeros(3)
                for iatom in range(self.natom):
                    R_center += np.sum(np.power(
                        Q_mat[3*iatom:3*(iatom+1), pmode], 2))*self.R_mat[iatom]
                dum += np.linalg.norm(R_center, ord=2)
        else:
            raise TypeError

        return dum

    def construct_unitary(self, pair_theta_dic):
        identity = np.eye(self.nmode)
        for pair, theta in pair_theta_dic.items():
            pmode = pair[0]
            qmode = pair[1]
            rotate_mat = np.eye(self.nmode)
            rotate_mat[pmode][pmode] = np.cos(theta)
            rotate_mat[qmode][qmode] = np.cos(theta)
            rotate_mat[pmode][qmode] = np.sin(theta)
            rotate_mat[qmode][pmode] = - np.sin(theta)
            identity = np.dot(identity, rotate_mat)
        return identity

    def run(self):
        def objective(theta, pmode, qmode):
            pair_theta = {}
            pair_theta[(pmode, qmode)] = theta
            rotate_matrix = self.construct_unitary(pair_theta)
            Q_mat = np.dot(self.Q_mat, rotate_matrix)
            return - self.metric(Q_mat, check_mode=[pmode, qmode])

        max_iter = 100
        pair_list = []
        for pmode, qmode in itertools.combinations(range(self.nmode), 2):
            if abs(self.freq[pmode] - self.freq[qmode]) < self.window:
                pair_list.append((pmode, qmode))
        pre_val = 0
        tol = 1.e-05
        met = self.metric(self.Q_mat)
        print('initial zeta = ', met)
        import random
        for k in range(max_iter):
            random.shuffle(pair_list)
            for pmode, qmode in pair_list:
                result = minimize_scalar(objective, bounds=(-np.pi, np.pi), 
                        args=(pmode, qmode))
                pair_theta = {}
                pair_theta[(pmode, qmode)] = result.x
                rotate_matrix = self.construct_unitary(pair_theta)
                self.Q_mat = np.dot(self.Q_mat, rotate_matrix)
                self.unitary = np.dot(self.unitary, rotate_matrix)
                met = self.metric(self.Q_mat)
            print(k, 'zeta = ', met, 'delta =', abs(pre_val - met))
            if abs(pre_val - met) < tol:
                break
            else:
                pre_val = met


        #def objective(array):
        #    pair_theta = {}
        #    cnt = 0
        #    #for pmode, qmode in itertools.combinations(range(self.nmode), 2):
        #    for pmode in range(self.nmode-1):
        #        #pair_theta[(pmode, qmode)] = array[cnt]
        #        pair_theta[(pmode, pmode+1)] = array[pmode]
        #        cnt += 1
        #    rotate_matrix = self.construct_unitary(pair_theta)
        #    Q_mat = np.dot(self.Q_mat, rotate_matrix)
        #    return - self.metric(Q_mat)

        #result = minimize(objective, np.zeros(self.nmode*(self.nmode-1)//2))
        #result = minimize(objective, np.zeros(self.nmode-1))
        #print(result.message)
        #print(result.fun)
        #if not result.success:
        #    assert False

        #pair_theta = {}
        #cnt = 0
        #for pmode, qmode in itertools.combinations(range(self.nmode), 2):
        #for pmode in range(self.nmode-1):
            #pair_theta[(pmode, qmode)] = result.x[cnt]
            #pair_theta[(pmode, pmode+1)] = result.x[pmode]
            #cnt += 1
        #rotate_matrix = self.construct_unitary(pair_theta)

        #self.Q_mat = np.dot(self.Q_mat, rotate_matrix)

        hess = np.dot(np.dot(
            self.unitary,  self.hess), self.unitary.T)

        return (self.Q_mat, list(np.sqrt(np.diag(hess))))

class GroupLocalizer:
    def __init__(self, vib, mw_hess: np.ndarray, domains: List[List[float]]):
        self.geom = vib.geom
        self.nmode = len(mw_hess)
        self.natom = len(self.geom)
        self.domains = domains
        self.mw_hess = mw_hess

    def run(self):
        self.unitary = np.zeros_like(self.mw_hess)
        for domain in self.domains:
            fancy_indices = []
            for iatom in domain:
                fancy_indices.extend([3*iatom, 3*iatom+1, 3*iatom+2])
            sub_hess = self.mw_hess[np.ix_(fancy_indices, fancy_indices)]
            _, v = eigh(sub_hess)
            self.unitary[np.ix_(fancy_indices, fancy_indices)] = v

        return (self.unitary, 
                np.sqrt(np.diag(self.unitary.T@self.mw_hess@self.unitary).tolist()))
