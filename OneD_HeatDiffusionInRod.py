import matplotlib.pyplot as plt
import numpy as np


class TempDistribution:
    def __init__(self, no_of_CV=5, len_of_rod=10.0, kw=400.0, ke=400.0, source=[3.0, 0.0, -4.0], \
                 T_x_w=373.0, T_x_e=573.0, heat_flux=10.0):

        self.no_of_CV = no_of_CV
        self.len_of_rod = len_of_rod
        self.kw = kw
        self.ke = ke
        # Linear source term of the form S = Sc+SpTp where Sp should be (-)ve
        self.Sc1 = source[0]
        self.Sc2 = source[1]
        self.Sp = source[2]
        # ....................
        self.T_x_w = T_x_w
        self.T_x_e = T_x_e
        self.q = heat_flux

    def temp_distribution(self, Tguess=None, BC='Drichlet'):

        delta_x = self.len_of_rod / self.no_of_CV  # width of each inner CV
        delta_x_e = delta_x
        delta_x_w = delta_x

        no_of_grid_points = self.no_of_CV + 2
        inner_grid_points = self.no_of_CV

        solution = np.zeros(inner_grid_points)

        noConvergence = True

        # assigning some default value to Tguess if it is not provided by the user
        if Tguess == None:
            Tguess = np.zeros(inner_grid_points)

        # checking the value of Sp for consistency
        if self.Sp <= 0.0:

            while (noConvergence == True):

                # these values and data structures will be re-initialised in every iteration
                matrix = []
                residue = np.zeros(inner_grid_points)
                b = np.zeros(inner_grid_points)
                convergence_count = 0
                # .......................................................

                for i in range(0, self.no_of_CV):  # loop over CVs

                    a = np.zeros(inner_grid_points)

                    if i == 0:  # west side outer CV

                        a[i + 1] = (self.ke / delta_x_e) * (-1)
                        a[i] = (self.ke / delta_x_e) + (self.kw / (delta_x_w / 2.0)) - self.Sp * delta_x

                        b[i] = (self.Sc1 + self.Sc2 * Tguess[i]) * delta_x + (self.kw / (delta_x_w / 2.0)) * self.T_x_w

                    elif i == (inner_grid_points - 1):  # east side outer CV

                        if BC == 'Drichlet':

                            a[i - 1] = (self.kw / delta_x_w) * (-1)
                            a[i] = (self.kw / delta_x_w) + (self.ke / (delta_x_e / 2.0)) - self.Sp * delta_x

                            b[i] = (self.Sc1 + self.Sc2 * Tguess[i]) * delta_x + (self.ke / (
                            delta_x_e / 2.0)) * self.T_x_e

                        elif BC == 'Neumann':

                            del_x = delta_x
                            del_x_w = delta_x_w

                            if self.q != 0.0:  # then divide the east side last CV into half to capture gradient in temp.
                                del_x = delta_x / 2.0
                                del_x_w = delta_x_w / 4.0

                            a[i - 1] = (self.kw / del_x_w) * (-1)
                            a[i] = (self.kw / del_x_w) - (self.Sp * del_x)

                            b[i] = (self.Sc1 + self.Sc2 * Tguess[i]) * (del_x) - self.q

                    else:  # inner CVs

                        a[i - 1] = (self.kw / delta_x_w) * (-1)
                        a[i] = (self.kw / delta_x_w) + (self.ke / delta_x_e) - self.Sp * delta_x
                        a[i + 1] = (self.ke / delta_x_e) * (-1)

                        b[i] = (self.Sc1 + self.Sc2 * Tguess[i]) * delta_x

                    matrix.append(a)

                solution = np.linalg.solve(np.array(matrix), b)  # NumPy method to solve system of linear eqns.

                # calculating residues in every iteration for each grid point within CVs
                for j in range(inner_grid_points):
                    residue[j] = abs(solution[j] - Tguess[j])
                # print (residue)

                # checking convergence.........
                for j in range(len(residue)):

                    if residue[j] < 0.01:
                        convergence_count += 1

                if convergence_count == inner_grid_points:

                    noConvergence = False

                else:

                    Tguess = solution  # now solution becomes Tguess for the next iteration
                    # ............................
        else:
            raise ValueError('Value of Sp should be negative for physical consistency of the discretisation')

        return (solution)








t_d = TempDistribution(no_of_CV=5, len_of_rod=1.0, kw=1.0, ke=1.0, source=[500.0, 0.0, -25.0], T_x_w=100.0,
                       heat_flux=0.0)
solution = t_d.temp_distribution(BC='Neumann')
# solution=solution[:]-273.0
print(solution)
X = np.linspace(0, 1, 5)
Y = solution
plt.plot(X, Y)
plt.show()
