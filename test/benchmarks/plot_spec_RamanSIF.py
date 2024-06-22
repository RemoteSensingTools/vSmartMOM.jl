import matplotlib.pyplot as plt
import numpy as np

# Step 2: Generate sample data
x = np.linspace(0, 10, 100)
y1 = np.sin(x)
y2 = np.cos(x)
y3 = np.tan(x)
y4 = np.exp(-x)

# Step 3: Create figure and subplots
fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(14, 10))

# Step 4: Plot data on each subplot
for i, ax in enumerate(axs.flat):
    ax.set_title(f"Plot {i+1}")
    if i < 4:
        ax.plot(x, [y1, y2, y3, y4][i], label=["sin(x)", "cos(x)", "tan(x)", "exp(-x)"][i])
    else:
        ax.plot(x + 5, [2*y1, 0.5*y2, 0.8*y3, 0.3*y4][i-4], label=["2*sin(x)", "0.5*cos(x)", "0.8*tan(x)", "0.3*exp(-x)"][i-4])
    ax.set_xlim(0, 15)
    ax.set_ylim(-1.5, 1.5)
    ax.legend()

# Step 5: Adjust layout and display
plt.tight_layout()
plt.show()
