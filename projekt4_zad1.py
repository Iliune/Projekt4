import numpy as np
import matplotlib.pyplot as plt

beta = 1.0
sigma = 1.0
gamma = 0.1

s_0, e_0, i_0, r_0 = 0.99, 0.01, 0, 0

def pochodne(s, e, i, r):
    ds = -beta * i * s
    de = beta * i * s - sigma * e
    di = sigma * e - gamma * i
    dr = gamma * i
    return ds, de, di, dr

def metoda_rungego_kutty(s, e, i, r, h):
    k1_s, k1_e, k1_i, k1_r = pochodne(s, e, i, r)
    k2_s, k2_e, k2_i, k2_r = pochodne(s + 0.5 * h * k1_s, e + 0.5 * h * k1_e, i + 0.5 * h * k1_i, r + 0.5 * h * k1_r)
    k3_s, k3_e, k3_i, k3_r = pochodne(s + 0.5 * h * k2_s, e + 0.5 * h * k2_e, i + 0.5 * h * k2_i, r + 0.5 * h * k2_r)
    k4_s, k4_e, k4_i, k4_r = pochodne(s + h * k3_s, e + h * k3_e, i + h * k3_i, r + h * k3_r)

    s_new = s + (h / 6) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s)
    e_new = e + (h / 6) * (k1_e + 2 * k2_e + 2 * k3_e + k4_e)
    i_new = i + (h / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i)
    r_new = r + (h / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)

    return s_new, e_new, i_new, r_new

t_max = 50
h = 0.1
steps = int(t_max / h) + 1

czas = np.linspace(0, t_max, steps)
podatni, narażeni, zainfekowani, martwi = np.zeros(steps), np.zeros(steps), np.zeros(steps), np.zeros(steps)

podatni[0], narażeni[0], zainfekowani[0], martwi[0] = s_0, e_0, i_0, r_0

for i in range(1, steps):
    podatni[i], narażeni[i], zainfekowani[i], martwi[i] = metoda_rungego_kutty(podatni[i - 1], narażeni[i - 1], zainfekowani[i - 1], martwi[i - 1], h)

plt.plot(czas, podatni, label="Podatni", color='lightpink')
plt.plot(czas, narażeni, label="Narażeni", color='indianred')
plt.plot(czas, zainfekowani, label="Zainfekowani", color='royalblue')
plt.plot(czas, martwi, label="Usunięci", color='mediumvioletred')
plt.xlabel("Czas")
plt.ylabel("Populacja")
plt.legend()
plt.show()
