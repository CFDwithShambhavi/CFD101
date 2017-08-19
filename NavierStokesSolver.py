import matplotlib.pyplot as plt
import numpy as np

class TwoDNavierStokes:

    def __init__(self):

        self.matrix = []

    def solver(self, no_of_cv, initial_p, initial_u, initial_v, px_index, py_index, py_val, px_val, u_index, u_val, v_index, v_val):

        analysis = '2D'
        length = 10.0  # in meter
        row = 2.0
        mu = 1.0
        g = 9.8
        delta_t = 1.0

        u_grid_dim = (5, 4)
        u_rows, u_cols = u_grid_dim
        v_grid_dim = (u_grid_dim[1], u_grid_dim[0])
        v_rows, v_vols = v_grid_dim
        p_grid_dim = (u_grid_dim[0], u_grid_dim[0])
        p_rows, p_cols = p_grid_dim

        no_of_cv = (u_rows + u_cols - 3) / 2
        delta_x = length / no_of_cv
        delta_y = length / no_of_cv

        rows = u_rows + u_cols
        cols = u_rows + u_cols

        matrix = np.zeros((rows, cols))

        print(matrix)

        # given B.Cs.
        inlet_u = ('given', 76.0)
        outlet_u = ('not_given', 11.0)
        inlet_v = ('not_given', 5.0)
        outlet_v = ('not_given', 0.0)
        inlet_p = ('given', 12.0)
        outlet_p = ('not_given', 2.0)

        # prescribing B.Cs. on grid layout
        # pressure
        matrix[0::2, 0] = inlet_p[1]  # at inlet
        matrix[0::2, cols - 1] = outlet_p[1]  # at outlet
        # u velocity
        matrix[0::2, 1] = inlet_u[1]  # inlet
        matrix[0::2, cols - 2] = outlet_u[1]  # outlet
        # v velocity
        matrix[1::2, 0] = inlet_v[1]  # inlet
        matrix[1::2, cols - 1] = outlet_v[1]  # outlet

        print(matrix)

        # u velocity coefficients
        Fe = np.zeros(u_grid_dim)
        Fw = np.zeros(u_grid_dim)
        Fn = np.zeros(u_grid_dim)
        Fs = np.zeros(u_grid_dim)
        De = np.zeros(u_grid_dim)
        Dw = np.zeros(u_grid_dim)
        Dn = np.zeros(u_grid_dim)
        Ds = np.zeros(u_grid_dim)
        aP_o = np.zeros(u_grid_dim)
        delta_p = np.zeros(u_grid_dim)

        zero = np.zeros(u_grid_dim)

        Fw[:, 1:] = row / 2.0 * (matrix[0::2, 3::2] + matrix[0::2, 1:cols - 2:2])
        Fe[:, 0:u_cols - 1] = row / 2.0 * (matrix[0::2, 1:cols - 2:2] + matrix[0::2, 3::2])
        Fn[1:, :] = row / 2.0 * (matrix[2::2, 1::2] + matrix[0:cols - 1:2, 1::2])
        Fs[0:u_rows - 1, :] = row / 2.0 * (matrix[0:cols - 1:2, 1::2] + matrix[2::2, 1::2])

        De[:, 0:u_cols - 1] = mu / delta_x
        Dw[:, 1:] = mu / delta_x
        Dn[1:, :] = mu / delta_y
        Ds[0:u_rows - 1, :] = mu / delta_y

        delta_p = (matrix[0::2, 0:cols - 1:2] - matrix[0::2, 2::2]) / delta_x
        b = delta_p + row * g

        aE = (De + np.fmax(-Fe, zero)) / delta_x
        aW = (Dw + np.fmax(Fw, zero)) / delta_x
        aN = (Dn + np.fmax(-Fn, zero)) / delta_y
        aS = (Ds + np.fmax(Fs, zero)) / delta_y
        aP_o[:, :] = row / delta_t

        aP = aP_o + np.fmax(Fe, zero) / delta_x + np.fmax(-Fw, zero) / delta_x + np.fmax(Fn, zero) / delta_y + np.fmax(
            -Fs, zero) / delta_y + De / delta_x + Dw / delta_x + Dn / delta_y + Ds / delta_y

        print(aE, '\n')
        print(aW, '\n')
        print(aN, '\n')
        print(aS, '\n')
        print(aP, '\n')
        print(b)

        # Re-arranging the equations in the form: -aEUE + aPUP - aSUS = u_source, where u_source = aEUE + aWUW + b
        # After re-arranging in the above form, equations are solved using TDMA scheme for every North-South direction in the grid,
        # sweeping from West to East.

        u_values = matrix[0::2, 1::2]

        while True:

            u_vals = u_values.copy()
            u_source = np.zeros(u_rows)

            for col in range(0, u_cols):

                tdma = np.diag(aP[:, col], k=0) + np.diag(-aN[1:, col], k=1) + np.diag(-aS[0:u_rows - 1, col], k=-1)
                if col == 0:
                    u_source[:] = b[:, col] + u_values[:, col + 1]
                elif col == u_cols - 1:
                    u_source[:] = b[:, col] + u_values[:, col - 1]
                else:
                    u_source[:] = b[:, col] + u_values[:, col + 1] + u_values[:, col - 1]

                # implementing u boundary conditions on u-grid
                if inlet_u[0] == 'given' and col == 0 or outlet_u[0] == 'given' and col == u_cols - 1:
                    tdma[:] = 0.0
                    I = np.ones(u_rows)
                    tdma += np.diag(I, k=0)
                    if col == 0:
                        u_source[:] = inlet_u[1]
                    else:
                        u_source[:] = outlet_u[1]

                u_solution = np.linalg.solve(tdma, u_source)
                u_values[:, col] = u_solution[:]

            residue = abs(u_vals[:] - u_values[:])
            print('************ Residue **************', '\n')
            print(residue, '\n')
            if residue[residue < 0.01].shape == (u_rows * u_cols,):
                break

        print('*************** u-velocity solution ****************')
        print(u_values, '\n')
        # print(u_solution)
        # print(tdma, ' = ', source)
