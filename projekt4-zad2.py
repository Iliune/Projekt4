import numpy as np
import matplotlib.pyplot as plt

beta, sigma, gamma = 0.5, 1.0, 0.1
s_0, e_0, i_0, r_0 = 0.99, 0.01, 0, 0

def pochodne(s, e, i, r):
    ds = -beta * i * s
    de = beta * i * s - sigma * e
    di = sigma * e - gamma * i
    dr = gamma * i
    return ds, de, di, dr

def runge_kutta(s, e, i, r, h):
    k1_s, k1_e, k1_i, k1_r = pochodne(s, e, i, r)
    k2_s, k2_e, k2_i, k2_r = pochodne(s + 0.5 * h * k1_s, e + 0.5 * h * k1_e, i + 0.5 * h * k1_i, r + 0.5 * h * k1_r)
    k3_s, k3_e, k3_i, k3_r = pochodne(s + 0.5 * h * k2_s, e + 0.5 * h * k2_e, i + 0.5 * h * k2_i, r + 0.5 * h * k2_r)
    k4_s, k4_e, k4_i, k4_r = pochodne(s + h * k3_s, e + h * k3_e, i + h * k3_i, r + h * k3_r)

    s_nowe = s + (h / 6) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s)
    e_nowe = e + (h / 6) * (k1_e + 2 * k2_e + 2 * k3_e + k4_e)
    i_nowe = i + (h / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i)
    r_nowe = r + (h / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)

    return s_nowe, e_nowe, i_nowe, r_nowe

t_max, h = 50, 0.1
liczba_krokow = int(t_max / h) + 1

czas = np.linspace(0, t_max, liczba_krokow)
podatni, narażeni, zainfekowani, martwi = np.zeros(liczba_krokow), np.zeros(liczba_krokow), np.zeros(liczba_krokow), np.zeros(liczba_krokow)

podatni[0], narażeni[0], zainfekowani[0], martwi[0] = s_0, e_0, i_0, r_0

for i in range(1, liczba_krokow):
    podatni[i], narażeni[i], zainfekowani[i], martwi[i] = runge_kutta(podatni[i - 1], narażeni[i - 1], zainfekowani[i - 1], martwi[i - 1], h)

plt.plot(czas, podatni, label="Podatni", color='lightpink')
plt.plot(czas, narażeni, label="Narażeni", color='indianred')
plt.plot(czas, zainfekowani, label="Zainfekowani", color='royalblue')
plt.plot(czas, martwi, label="Martwi", color='mediumvioletred')
plt.xlabel("Czas")
plt.ylabel("Populacja")
plt.legend()
plt.show()

'''
Wartość współczynnika beta wpływa na tempo rozprzestrzeniania się wirusa w populacji. 
Niższa wartość beta skutkuje wolniejszym tempem rozprzestrzeniania, mniejszą liczbą zakażeń jednocześnie, 
ale również dłuższą trwająca epidemią. 
Wyższa wartość beta oznacza szybsze rozprzestrzenianie się wirusa,
większą liczbę zarażonych w krótszym czasie, potencjalnie większe obciążenie systemu zdrowia, 
ale także krótszą trwającą epidemią. 
Wybór wartości beta wpływa na balans między kontrolą epidemii a minimalizacją skutków ubocznych.
'''
