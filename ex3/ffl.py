import numpy as np
from plotly import graph_objects as go

Y_B_MAX = 10
X_Y_KD = 0.8
X_Z_KD = 0.7
Y_ALPHA = 0.4
Z_B_MAX = 10
Y_Z_KD = 0.8
Z_ALPHA = 0.6
X_ALPHA = 0.6
X_STST = 8


def dy_a(x, y):
    return Y_B_MAX * (x / (x + X_Y_KD)) - y * Y_ALPHA


def dz_a(x, y, z):
    numinator = (x / X_Z_KD) + (y / Y_Z_KD)
    return Z_B_MAX * (numinator / (1 + numinator)) - (z * Z_ALPHA)

def dz_c(x, y, z):
    return Z_B_MAX * ((x / (x + X_Z_KD)) * (y / (y + Y_Z_KD))) - z * Z_ALPHA


def dz_d(x, y, z):
    return Z_B_MAX * ((x / (x + X_Z_KD)) - (x / (x + X_Z_KD)) * (y / (y + Y_Z_KD))) - z * Z_ALPHA


def iteration_a(x, y, z, dt, signal):
    if not signal:
        x = 0
    y_new = y + dy_a(x, y) * dt
    z_new = z + dz_a(x, y, z) * dt
    return y_new, z_new


def dx_b(prev_x, signal, t):
    if signal:
        return X_STST * (1 - np.exp(-X_ALPHA * t))
    else:
        return min(prev_x, X_STST * np.exp(-X_ALPHA * t))

def iteration_b(x, y, z, dt, signal, t):
    x_new = dx_b(x, signal, t)
    if not signal:
        x = 0
    y_new = y + dy_a(x, y) * dt
    z_new = z + dz_a(x, y, z) * dt
    return x_new, y_new, z_new


def iteration_c(x, y, z, dt, signal):
    if not signal:
        x = 0
    y_new = y + dy_a(x, y) * dt
    z_new = z + dz_c(x, y, z) * dt
    return y_new, z_new


def iteration_d(x, y, z, dt, signal):
    if not signal:
        x = 0
    y_new = y + dy_a(x, y) * dt
    z_new = z + dz_d(x, y, z) * dt
    return y_new, z_new



def plot_system(x, y, z, name, dt=0.1):
    # plot the system
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=np.arange(0, 100, dt), y=y, name="Y"))
    fig.add_trace(go.Scatter(x=np.arange(0, 100, dt), y=z, name="Z"))
    fig.add_trace(go.Scatter(x=np.arange(0, 100, dt), y=x, name="X"))
    # add signal line
    fig.add_shape(type="line",
                  x0=20*dt, y0=0, x1=20*dt, y1=15,
                  line=dict(color="black", width=2, dash="dash"))
    fig.update_layout(title=name, xaxis_title="Time", yaxis_title="Concentration")
    fig.show()
    file_name = f"ffl_{name}.png"
    fig.write_image(format="png", engine="kaleido", file=file_name)


if __name__ == '__main__':
    # A
    x = np.ones(100) * X_STST
    y = np.zeros(100)
    z = np.zeros(100)
    dt = 0.1
    for i in range(1, 100):
        if i <= 20:
            signal = True
        else:
            signal = False
        y[i], z[i] = iteration_a(x[i], y[i - 1], z[i - 1], dt, signal)

    plot_system(x, y, z, "A", dt)

    # B
    x = np.zeros(100)
    y = np.zeros(100)
    z = np.zeros(100)
    dt = 0.1

    signal_end = 20
    for i in range(1, 100):
        if i <= signal_end:
            signal = True
            t = i * dt
        else:
            signal = False
            t = (i - signal_end) * dt

        x[i], y[i], z[i] = iteration_b(x[i - 1], y[i - 1], z[i - 1], dt, signal, t)

    plot_system(x, y, z, "B", dt)

    # C
    x = np.ones(100) * X_STST
    y = np.zeros(100)
    z = np.zeros(100)
    dt = 0.1

    for i in range(1, 100):
        if i <= 20:
            signal = True
        else:
            signal = False
        y[i], z[i] = iteration_c(x[i - 1], y[i - 1], z[i - 1], dt, signal)

    plot_system(x, y, z, "C", dt)

    # D
    x = np.ones(100) * X_STST
    y = np.zeros(100)
    z = np.zeros(100)
    dt = 0.1

    for i in range(1, 100):
        if i <= 20:
            signal = True
        else:
            signal = False
        y[i], z[i] = iteration_d(x[i - 1], y[i - 1], z[i - 1], dt, signal)

    plot_system(x, y, z, "D", dt)







