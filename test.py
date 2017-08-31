import numpy as np

def display(var = ((0.0,0.0),)):
    #print (var)
    for i in var:
        if i[0] != 0.0:
            a,b = i
            #print(a+b, '\n')

    no_of_cv = 5
    no_of_grid_points = no_of_cv * 2 + 1
    p_u_grid = np.zeros((no_of_grid_points, no_of_grid_points))
    initial_p = 100.0
    initial_u = 18.0
    # initial values to velocity and pressure
    p_u_grid[1::2][1::2] = initial_p
    p_u_grid[0::2][0::2] = initial_u

    print(p_u_grid)

display(((2.0,1.0),(5.0,6.0),(8.0,7.0)))