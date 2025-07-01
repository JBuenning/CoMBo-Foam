import numpy as np
import matplotlib.pyplot as plt
import os


def get_time_directories(path):
    times = []
    for directory in os.listdir(path):
        try:
            float(directory)
            times.append(directory)
        except: pass
    times.sort(key=float)
    return times

def extract_from_time(time_name):
    point_file = os.path.join(time_name, "polyMesh/points")
    points = []
    with open(point_file) as f:
        for l in f.readlines():
            if len(l) > 3 and l[0] == '(':
                candidate = l.strip('\n()')
                x, y, z = candidate.split()
                points.append([float(x), float(y), float(z)])

    tol = 1e-5
    points = np.array(points)
    down_mask = points[:,2] < 0.5
    right_mask = down_mask & (points[:,0] > 1 - tol) & (points[:,0] < 1 + tol)
    left_mask = down_mask & (points[:,0] > 0 - tol) & (points[:,0] < 0 + tol)

    right_value = np.max(points[right_mask][:,1])
    left_value = np.max(points[left_mask][:,1])

    return left_value, right_value


def extract_points(path):
    times_str = sorted(get_time_directories(path))
    times_numeric = [float(t) for t in times_str]
    assert times_str[0] == '0'
    times_str = ["constant", *times_str[1:]]
    left = []
    right = []

    for t_str in times_str:
        l, r = extract_from_time(os.path.join(path, t_str))
        left.append(l)
        right.append(r)

    return np.array(times_numeric), np.array(left), np.array(right)

t_other = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
left_gupta_A = [0.4144, 0.5016, 0.5756, 0.6410, 0.7002, 0.7547, 0.8055, 0.8532]
left_gupta_B = [0.4084, 0.4938, 0.5664, 0.6308, 0.6892, 0.7429, 0.7928, 0.8396]
right_gupta_A = [0.5533, 0.6085, 0.6623, 0.7135, 0.7621, 0.8083, 0.8524, 0.8945]
right_gupta_B = [0.5460, 0.5990, 0.6509, 0.7003, 0.7472, 0.7918, 0.8344, 0.8752]

fig, ax = plt.subplots(layout='constrained',)

t, left, right = extract_points(path='.')
ax.plot(t, left, color='tab:blue', label=r'Simulation, $x_1 = 0$')
ax.scatter(t_other, left_gupta_A, color='tab:blue', marker='x', label=r'Gupta (A), $x_1 = 0$')
ax.scatter(t_other, left_gupta_B, color='tab:blue', marker='^', label=r'Gupta (B), $x_1 = 0$')

ax.plot(t, right, color='tab:orange', label=r'Simulation, $x_1 = 1$')
ax.scatter(t_other, right_gupta_A, color='tab:orange', marker='x', label=r'Gupta (A), $x_1 = 1$')
ax.scatter(t_other, right_gupta_B, color='tab:orange', marker='^', label=r'Gupta (B), $x_1 = 1$')
ax.grid()
ax.set_xlim(left=0, right=0.42)
ax.set_ylim(top=0.935)
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$s(x_1, t)$")
plt.legend()
plt.show()

