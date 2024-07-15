import uproot
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import library as lib

filename = "scanx.200mrad.A1.txt.tree.root"

def file_to_list(filename):
    file = uproot.open(filename)
    tree = file["tree"]
    data = tree.arrays(["x", "y", "rate"], library="np")
    x = data["x"]
    y = data["y"]
    rate = data["rate"]
    return x, y, rate

# Open the file and read the data
x, y, rate = file_to_list(filename)

# Remove zero values from rate and corresponding x values
nonzero_mask = rate != 0
x_nonzero = x[nonzero_mask]
rate_nonzero = rate[nonzero_mask]

print(f"Original length of x: {len(x)}, Nonzero length of x: {len(x_nonzero)}")
print(f"Original length of rate: {len(rate)}, Nonzero length of rate: {len(rate_nonzero)}")

# Calculate average rate and error for each unique x value
rate_dict = defaultdict(list)
for x_val, rate_val in zip(x_nonzero, rate_nonzero):
    rate_dict[x_val].append(rate_val)

avg_x = []
avg_rate = []
rate_errors = []

for x_val in sorted(rate_dict.keys()):
    rates = rate_dict[x_val]
    avg_x.append(x_val)
    avg_rate.append(np.mean(rates))
    rate_errors.append(np.std(rates) / np.sqrt(len(rates)))


print(len(avg_x), len(avg_rate))
# Plotting
dx, drate = lib.derivative(avg_x, avg_rate)
plt.errorbar(avg_x, avg_rate, yerr=rate_errors, fmt='r.', label="Average Rate [Hz]")
plt.plot( dx, drate, label = " derivative")
plt.xlabel("x")
plt.ylabel("Average Rate [Hz]")
plt.title("Average Rate vs x with Error Bars")
plt.legend()
plt.show()
