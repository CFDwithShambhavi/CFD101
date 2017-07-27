import matplotlib.pyplot as plt
import numpy as np

class TempDistribution:

    def __init__(self, no_of_CV, len_of_rod, kw, ke, source, T_x_w, T_x_e):

        self.no_of_CV = no_of_CV
        self.len_of_rod = len_of_rod
        self.kw = kw
        self.ke = ke
        #Linear source term of the form S = Sc + SpTp
        self.Sc = source[0]
        self.Sp = source[1]
        #....................
        self.T_x_w = T_x_w
        self.T_x_e = T_x_e

    def temp_distribution(self):

        delta_x = self.len_of_rod / self.no_of_CV #width of each inner CV
        delta_x_e = delta_x
        delta_x_w = delta_x
        no_of_grid_points = self.no_of_CV + 2
        inner_grid_points = self.no_of_CV
        b = np.zeros(inner_grid_points)
        matrix = []

        for i in range(0,self.no_of_CV): #loop over CVs

            a = np.zeros(inner_grid_points)

            if i == 0: #west side outer CV

                a[i+1] = (self.ke / delta_x_e) * (-1)
                a[i] = (self.ke / delta_x_e) + (self.kw / (delta_x_w / 2.0)) + self.Sp * delta_x

                b[i] = self.Sc * delta_x + (self.kw / (delta_x_w / 2.0)) * self.T_x_w

            elif i == (inner_grid_points-1): #east side outer CV

                a[i-1] = (self.kw / delta_x_w) * (-1)
                a[i] = (self.kw / delta_x_w) + (self.ke / (delta_x_e / 2.0)) + self.Sp * delta_x

                b[i] = self.Sc * delta_x + (self.ke / (delta_x_e / 2.0)) * self.T_x_e

            else: #inner CVs

                a[i-1] = (self.kw / delta_x_w) * (-1)
                a[i] = (self.kw / delta_x_w) + (self.ke / delta_x_e) + self.Sp * delta_x
                a[i+1] = (self.ke / delta_x_e) * (-1)

                b[i] = self.Sc * delta_x

            matrix.append(a)

        solution = np.linalg.solve(np.array(matrix),b)

        return(solution)




t_d  = TempDistribution(5,0.02,0.5,0.5,[1000.0e+03,0.0],373.0,473.0)
solution = t_d.temp_distribution()
solution = solution[:] - 273.0
print(solution)
X = np.linspace(0,2,5)
Y = solution
plt.plot(X,Y)
plt.show()
