import matplotlib.pyplot as plt
import numpy as np


class TempDistribution:
    def __init__(self, no_of_CV=0, len_of_rod=0.0, kw=0.0, ke=0.0, source=(0.0, 0.0, 0.0),
                 T=((0.0, 0),), q=((0.0, 0),), h=0.0, amb_T=0.0, row=0.0, cp=0.0):

        self.no_of_CV = no_of_CV
        self.len_of_rod = len_of_rod
        self.kw = kw
        self.ke = ke
        self.h = h  # heat transfer coeff.
        self.row = row
        self.cp = cp
        self.amb_T = amb_T  # ambient temp.
        # Linear source term of the form S = Sc+SpTp where Sp should be (-)ve
        self.Sc1, self.Sc2, self.Sp = source
        # ....................
        self.T = T
        self.q = q
        # ....................

    def temp_distribution(self, temp_t0=0.0, guess_temp=0.0, large_num=1e+10, time=1.0, delta_t=1.0):

        result = []

        delta_x = self.len_of_rod / self.no_of_CV  # width of each inner CV
        delta_x_e = delta_x
        delta_x_w = delta_x

        no_of_grid_points = self.no_of_CV + 2

        initial_T = np.zeros(no_of_grid_points)
        guess_T = np.zeros(no_of_grid_points)

        t = 0.0  # iterating variable for time loop

        # calculating the value of a_p^o: coeff. of time neighbour for every grid point
        a_po = self.row * self.cp * delta_x / delta_t

        # assigning initial value to the temperature array
        if temp_t0 != 0.0:
            initial_T += temp_t0

        # guess temperature at every grid point in case of temp. dependent heat source
        if guess_temp != 0.0:
            guess_T += guess_temp

        # checking the value of Sp for consistency
        if self.Sp <= 0.0:

            while t < time:

                print('!!!!!!!!!!!!  Time = ', t, ' sec. !!!!!!!!!!!!! ', '\n')

                while True:

                    # these values and data structures will be re-initialised in every iteration
                    matrix = np.zeros((no_of_grid_points, no_of_grid_points))
                    residue = np.zeros(no_of_grid_points)
                    b = np.zeros(no_of_grid_points)
                    # .......................................................

                    for i in range(0, no_of_grid_points):  # loop over grid points within each CV

                        a = np.zeros(no_of_grid_points)

                        if i == 0:  # west side boundary grid point

                            a[i] = (self.ke / (delta_x_e / 2.0)) - self.Sp * (delta_x / 2.0) + a_po
                            a[i + 1] = (self.ke / (delta_x_e / 2.0)) * (-1)

                            b[i] = (self.Sc1 + self.Sc2 * guess_T[i]) * (delta_x / 2.0) + a_po * initial_T[i]

                        elif i == (no_of_grid_points - 1):  # east side boundary grid point

                            a[i] = (self.kw / (delta_x_w / 2.0)) - self.Sp * (delta_x / 2.0) + a_po
                            a[i - 1] = (self.kw / (delta_x_w / 2.0)) * (-1)

                            b[i] = (self.Sc1 + self.Sc2 * guess_T[i]) * (delta_x / 2.0) + a_po * initial_T[i]

                        else:  # inner grid points within inner CVs

                            del_x_w = delta_x_w
                            del_x_e = delta_x_e

                            if i == 1:
                                del_x_w = delta_x_w / 2.0
                            elif i == no_of_grid_points - 2:
                                del_x_e = delta_x_e / 2.0

                            a[i - 1] = (self.kw / del_x_w) * (-1)
                            a[i] = (self.kw / del_x_w) + (self.ke / del_x_e) - self.Sp * delta_x + a_po
                            a[i + 1] = (self.ke / del_x_e) * (-1)

                            b[i] = (self.Sc1 + self.Sc2 * guess_T[i]) * delta_x + a_po * initial_T[i]

                        matrix[i] = a[:]

                    # boundary condition implementation:
                    for var in self.T:
                        temp_val, grid_val = var
                        if temp_val != 0.0: # temperatures are specified at the boundaries
                            matrix[grid_val][grid_val] += large_num
                            b[grid_val] = temp_val * (large_num / matrix[grid_val][grid_val])
                            matrix[grid_val] /= matrix[grid_val][grid_val]

                    for var in self.q:
                        q_val, grid_val = var
                        if grid_val == 0 or grid_val == (no_of_grid_points - 1):
                            if q_val != 0.0: # direct implementation of heat flux at boundaries
                                b[grid_val] = b[grid_val] + q_val
                            elif q_val == 0.0 and self.h != 0.0 and self.amb_T != 0.0:# heat flux given as a fn of heat transfer coeff.
                                b[grid_val] = b[grid_val] + self.h * self.amb_T
                                matrix[grid_val][grid_val] += self.h
                        else:
                            raise ValueError('heat flux can only be applied at outermost grid points !!')

                    solution = np.linalg.solve(matrix, b)  # NumPy method to solve system of linear eqns.

                    # calculating residues in every iteration for each grid point within CVs
                    residue[:] = abs(solution[:] - guess_T[:])

                    print('Residue : ', residue)

                    # checking convergence.........
                    if all(residue < 0.01):

                        print('!!!!!!!!!!!! Convergence !!!!!!!!!!!!! ', '\n')
                        result.append(solution[:])
                        break

                    else:

                        guess_T = solution[:]  # now solution becomes guess_T for the next iteration
                        # ............................

                t += delta_t
        else:
            raise ValueError('Value of Sp should be negative for physical consistency of the discretisation !!')

        return result


# post processing
t_d = TempDistribution(no_of_CV=5, len_of_rod=0.02, kw=0.5, ke=0.5, source=(1000.0e+03, 0.0, 0.0),
                       T=((100.0, 0), (200.0, 6)), row=8960.0, cp=385.0)
sol = t_d.temp_distribution(time=5.0, delta_t=1.0)

X = np.linspace(0, 2.0, 7)

for i in range(0, len(sol)):
    plt.figure(1)
    print('!!!!!!!!!!!!!! solution !!!!!!!!!!!!', '\n')
    print(sol[i], '\n')
    Y = sol[i]
    plt.plot(X, Y)
    plt.show()
