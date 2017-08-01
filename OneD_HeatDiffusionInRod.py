import matplotlib.pyplot as plt
import numpy as np

class TempDistribution:
    def __init__(self, no_of_CV = 0, len_of_rod = 0.0, kw = 0.0, ke = 0.0, source = [0.0, 0.0, 0.0], \
                 bc = dict(T1 = (0.0,0), T2 = (0.0,0), q = (0.0)), h = 0.0, amb_T = 0.0):

        self.no_of_CV = no_of_CV
        self.len_of_rod = len_of_rod
        self.kw = kw
        self.ke = ke
        self.h = h #heat transfer coeff.
        self.amb_T = amb_T #ambient temp.
        # Linear source term of the form S = Sc+SpTp where Sp should be (-)ve
        self.Sc1 = source[0]
        self.Sc2 = source[1]
        self.Sp = source[2]
        # ....................
        self.bc = bc

    def temp_distribution(self, guess_T = None, large_num = 1e+10):

        global matrix
        delta_x = self.len_of_rod / self.no_of_CV  # width of each inner CV
        delta_x_e = delta_x
        delta_x_w = delta_x

        inner_grid_points = self.no_of_CV
        no_of_grid_points = self.no_of_CV + 2

        solution = np.zeros(no_of_grid_points)

        # assigning some default value to Tguess if it is not provided by the user
        if guess_T == None:
            guess_T = np.zeros(no_of_grid_points)

        # checking the value of Sp for consistency
        if self.Sp <= 0.0:

            while (True):

                # these values and data structures will be re-initialised in every iteration
                matrix = np.zeros((no_of_grid_points, no_of_grid_points))
                residue = np.zeros(no_of_grid_points)
                b = np.zeros(no_of_grid_points)
                convergence_count = 0
                # .......................................................

                for i in range(0, no_of_grid_points):  # loop over grid points within each CV

                    a = np.zeros(no_of_grid_points)

                    if i == 0: # west side boundary grid point

                        a[i] = (self.ke / (delta_x_e / 2.0)) - self.Sp * (delta_x / 2.0)
                        a[i + 1] = (self.ke / (delta_x_e / 2.0)) * (-1)

                        b[i] = (self.Sc1 + self.Sc2 * guess_T[i]) * (delta_x / 2.0)

                    elif i == (no_of_grid_points - 1):  # east side boundary grid point

                        a[i] = (self.kw / (delta_x_w / 2.0)) - self.Sp * (delta_x / 2.0)
                        a[i - 1] = (self.kw / (delta_x_w / 2.0)) * (-1)

                        b[i] = (self.Sc1 + self.Sc2 * guess_T[i]) * (delta_x / 2.0)

                    else:  # inner grid points within inner CVs

                        del_x_w = delta_x_w
                        del_x_e = delta_x_e

                        if i == 1:
                            del_x_w = delta_x_w / 2.0
                        elif i == no_of_grid_points - 2:
                            del_x_e = delta_x_e / 2.0

                        a[i - 1] = (self.kw / del_x_w) * (-1)
                        a[i] = (self.kw / del_x_w) + (self.ke / del_x_e) - self.Sp * delta_x
                        a[i + 1] = (self.ke / del_x_e) * (-1)

                        b[i] = (self.Sc1 + self.Sc2 * guess_T[i]) * delta_x

                    matrix[i] = a[:]

                #boundary condition implementation:
                for k, v in self.bc.items():
                    if (k == 'T1' or k == 'T2') and v[0] != 0.0: # Drichlet
                        index = v[1]
                        matrix[index][index] += large_num
                        b[index] = v[0] * (large_num / matrix[index][index])
                        matrix[index] /= matrix[index][index]

                    elif k == 'q' and (self.bc['T1'][0] == 0.0 or self.bc['T2'][0] == 0.0):
                        index = v[1]
                        if index == 0 or index == (no_of_grid_points - 1): # Neumann
                            if v[0] != 0.0 and self.h == 0.0:
                                b[index] = b[index] + self.q

                            elif v[0] == 0.0 and self.h != 0.0 and self.amb_T != 0.0: # Mixed type
                                b[index] = b[index] + self.h * self.amb_T
                                matrix[index][index] += self.h

                        else:
                            raise ValueError('heat flux can only be applied at outermost grid points !!')

                solution = np.linalg.solve(matrix, b)  # NumPy method to solve system of linear eqns.

                # calculating residues in every iteration for each grid point within CVs
                residue[:] = abs(solution[:] - guess_T[:])

                print ('Residue : ',residue)

                # checking convergence.........
                if all(residue[:] < 0.01):

                    print ('Convergence !!!!! ')
                    break

                else:

                    guess_T = solution[:]  # now solution becomes guess_T for the next iteration
                    # ............................
        else:
            raise ValueError('Value of Sp should be negative for physical consistency of the discretisation !!')

        return (solution)







#post processing
t_d = TempDistribution(no_of_CV = 17, len_of_rod = 1.0, kw = 1.0, ke = 1.0, source = [500.0, 25.0, -50.0], \
                 bc = dict(T1 = (100.0,0), T2 = (0.0,16), q = (0.0,18)))
solution = t_d.temp_distribution()
print(solution)
X = np.linspace(0, 1.0, 19)
Y = solution[:]
plt.plot(X, Y)
plt.show()
