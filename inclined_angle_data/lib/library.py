import numpy as np
import matplotlib.pyplot as plt

def derivative(x, y):
    dydx = []
    x_axis = []
    if len(x) == len(y):
        for i in range(len(x) - 1):
            dy = y[i + 1] - y[i]
            dx = x[i + 1] - x[i]
            newX = (x[i + 1] + x[i]) / 2
            dydx.append(dy / dx)
            x_axis.append(newX)
    else:
        print("The size of the vectors are not the same")
    return np.array(x_axis), np.array(dydx)

