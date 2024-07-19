import numpy as np
import matplotlib.pyplot as plt

def circle_overlap(x, R, c, bg):
    return bg + R * R * np.arccos(1 - (x - c) / R) - (R - (x - c)) * np.sqrt(2 * R * (x - c) - (x - c) * (x - c))

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


filename = "scanx.200mrad.A1.txt"
data = np.loadtxt(filename, delimiter=' ')  # adjust delimiter if necessary

# Extract the 1st and 3rd columns
column1 = data[:, 0]  # first column
column3 = data[:, 2]  # third column

x_data = np.array(column1)
rate = np.array(column3)

R = 10
x = np.linspace(0, 2 * R, 1000)
#y = circle_overlap(x, R)


#calculating the derivatives
dydx, x_mid = derivative(x, y) #1st derivative
d2ydx2, x_mid_2 = derivative(x_mid, dydx) #2nd derivative


x_data1, drate = derivative(x_data, rate)
x_data2, ddrate = derivative(x_data1, drate)

plt.plot(column1, column3, "r.-", label = 'data')
plt.plot(x_data1, drate , "b.-", label = 'first Derivative')
#plt.plot(x_data2, ddrate , "g.-", label='Second Derivative')

# plt.plot(x, circle_overlap(x, R), "r.", label='circle_overlap(x)')
# plt.plot(x_mid, dydx, "*--", label='first Derivative')
# plt.plot(x_mid_2, d2ydx2, "+--", label='Second Derivative')
plt.xlabel("x")
plt.ylabel("y and derivatives")
plt.legend()
plt.show()
