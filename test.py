def display(var = ((0.0,0.0),)):
    print (var)
    for i in var:
        if i[0] != 0.0:
            a,b = i
            print(a+b, '\n')

display(((2.0,1.0),(5.0,6.0),(8.0,7.0)))