import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

# Step 1: Get the number of points
num_points = int(input("How many points? "))

# Step 2: Get the points from the user
points = []
print("Enter the points: (x y)")
for i in range(num_points):
    x_str, y_str = input().split()
    x, y = float(x_str), float(y_str)
    points.append((x, y))

# Convert points to numpy array
points = np.array(points)

# Step 3: Compute all pairwise distances
# Create a lookup table for distances
num_points = len(points)
distances = {}  # key: (i, j), value: distance between point i and point j

for i in range(num_points):
    for j in range(num_points):
        p_i = points[i]
        p_j = points[j]
        dist = np.linalg.norm(p_i - p_j)
        distances[(i, j)] = dist

# Step 4: Ask the user for starting point v0
print("Enter starting point: (x y)")
v0_str = input()
v0_x_str, v0_y_str = v0_str.split()
v0_x, v0_y = float(v0_x_str), float(v0_y_str)
v0 = np.array([v0_x, v0_y])

# Check that v0 is in points
v0_index = None
for idx, p in enumerate(points):
    if np.allclose(p, v0):
        v0_index = idx
        break

if v0_index is None:
    print("v0 is not in the set of points.")
    exit()

# Step 5: Compute S
S = set()
for idx, p in enumerate(points):
    if idx == v0_index:
        continue
    dist = distances[(v0_index, idx)]
    if dist == 0:
        continue  # Avoid log2(0)
    s_element = int(np.floor(-np.log2(dist)))
    for delta in range(-2, 3):  # -2, -1, 0, 1, 2
        S.add(s_element + delta)
S = sorted(S)
print(f"S = {S}")

# Step 6: Generate the nets
nets = []

# Initialize the first net with v0
net = [v0_index]
nets.append(net.copy())

# For s in S, generate the nets
for s in S:
    previous_net = net.copy()
    new_net = previous_net.copy()
    for idx in range(num_points):
        if idx in previous_net:
            continue
        # Check that the point is at least 2^{-s} away from all points in previous net
        val = True
        for p_idx in previous_net:
            dist = distances[(p_idx, idx)]
            if dist < 2 ** (-s):
                val = False
                break
        if val:
            new_net.append(idx)
    net = new_net
    # Check if the net is unique (not already in nets)
    if sorted(net) not in [sorted(n) for n in nets]:
        nets.append(net.copy())

# Now we have the list of nets in 'nets'
# Each net is a list of indices into 'points'

# Step 7: Create point_appearances
point_appearances = {}  # point index -> net_index where it first appears

for net_idx, net in enumerate(nets):
    for idx in net:
        if idx not in point_appearances:
            point_appearances[idx] = net_idx

# Step 8: Visualize the points and the nets
# We can use matplotlib to plot the points

# Create a figure and axes
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)  # Make room for the buttons

# Plot all the points in black initially
x_vals = points[:, 0]
y_vals = points[:, 1]
colors = ['black'] * num_points
sc = ax.scatter(x_vals, y_vals, c=colors)

# net_index indicates the current net being displayed
net_index = 0  # Start with the first net

def update_plot():
    # Update colors based on net_index
    colors = ['black'] * num_points
    for idx in range(num_points):
        if idx in point_appearances and point_appearances[idx] <= net_index:
            colors[idx] = 'C{}'.format(point_appearances[idx] % 10)
    sc.set_color(colors)
    ax.set_title(f"Net {net_index+1}/{len(nets)}")
    fig.canvas.draw_idle()

def next_net(event):
    global net_index
    if net_index < len(nets) - 1:
        net_index += 1
        update_plot()

def prev_net(event):
    global net_index
    if net_index > 0:
        net_index -= 1
        update_plot()

# Initial plot
update_plot()

# Add buttons to step through the nets
axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
btn_prev = Button(axprev, 'Previous')
btn_prev.on_clicked(prev_net)

axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
btn_next = Button(axnext, 'Next')
btn_next.on_clicked(next_net)

plt.show()
