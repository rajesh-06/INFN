import numpy as np
import matplotlib.pyplot as plt

def circle_segment_area(R, x):
    return R * R * np.arccos(1 - x / R) - (R - x) * np.sqrt(2 * R * x - x * x)

    #return R**2 * np.arccos(d/R) - d * np.sqrt(R**2 - d**2)

def overlap_area(a, R, x):
    if x < -R or x > a + R:
        return 0  # No overlap
    elif R <= x <= a - R:
        return np.pi * R**2  # Full overlap
    elif x > -R and x < R:
        return circle_segment_area(R, x)
    elif x > a-R and x < a+R:
        return circle_segment_area(R, a + R-x)
    
def circ(x, a, b, R):
    return b + np.sqrt(R*R - (x-a)**2)
#   #  elif -R <= x < R:
#  #       d = R - x
# #        return circle_segment_area(R, d) + (x + R) * a
#     elif a - R < x <= a:
#         d = x - (a - R)
#         return np.pi * R**2 - circle_segment_area(R, d) + (a - x + R) * a
#     else:
#         if x < a/2:
#             d1 = R - x
#             d2 = x + R - a
#         else:
#             d1 = x - (a - R)
#             d2 = a + R - x
#         return circle_segment_area(R, d1) + circle_segment_area(R, d2)

# Parameters
a = 20  # Side length of the square
R = 5   # Radius of the circle
x = np.linspace(-R, a + R, 1000)

# Calculate overlap area
# overlap_areas = [overlap_area(a, R, xi) for xi in x]

# # Plotting
# # plt.plot(x, overlap_areas, "r--")
# plt.plot(x, circ(x, 0, 2, R), "r--")

# plt.xlabel('Horizontal Distance (x)')
# plt.ylabel('Overlap Area')
# plt.title('Overlap Area of Circle Moving Over Square')
# plt.grid(True)
# plt.show()

list1 = [1,2, 0, 0, 1, 2, 345]
zero_ind = []
for i in range(len(list1)):
    if list1[i] == 0:
        zero_ind.append(i)
print(zero_ind)
for i in range(len(zero_ind),0,-1):
    print(i-1)
    print(zero_ind[i-1])
    del list1[zero_ind[i-1]]
print(list1)

