import numpy as np
import matplotlib.pyplot as plt

s_0 = 0.99
e_0 = 0.01
i_0 = 0
r_0 = 0
sigma = 1.0
gamma = 0.1

def podatni(s, e, i, beta):
    return -beta * i * s

def narażeni(s, e, i, beta, sigma):
    return beta * i * s - sigma * e

def zainfekowani(e, i, sigma, gamma):
    return sigma * e - gamma * i

def martwi(i, gamma):
    return gamma * i

def metoda_rungego_kutty(s, e, i, r, beta, sigma, gamma, h):
    k1_s = podatni(s, e, i, beta)
    k1_e = narażeni(s, e, i, beta, sigma)
    k1_i = zainfekowani(e, i, sigma, gamma)
    k1_r = martwi(i, gamma)

    k2_s = podatni(s + 0.5 * h * k1_s, e + 0.5 * h * k1_e, i + 0.5 * h * k1_i, beta)
    k2_e = narażeni(s + 0.5 * h * k1_s, e + 0.5 * h * k1_e, i + 0.5 * h * k1_i, beta, sigma)
    k2_i = zainfekowani(e + 0.5 * h * k1_e, i + 0.5 * h * k1_i, sigma, gamma)
    k2_r = martwi(i + 0.5 * h * k1_i, gamma)

    k3_s = podatni(s + 0.5 * h * k2_s, e + 0.5 * h * k2_e, i + 0.5 * h * k2_i, beta)
    k3_e = narażeni(s + 0.5 * h * k2_s, e + 0.5 * h * k2_e, i + 0.5 * h * k2_i, beta, sigma)
    k3_i = zainfekowani(e + 0.5 * h * k2_e, i + 0.5 * h * k2_i, sigma, gamma)
    k3_r = martwi(i + 0.5 * h * k2_i, gamma)

    k4_s = podatni(s + h * k3_s, e + h * k3_e, i + h * k3_i, beta)
    k4_e = narażeni(s + h * k3_s, e + h * k3_e, i + h * k3_i, beta, sigma)
    k4_i = zainfekowani(e + h * k3_e, i + h * k3_i, sigma, gamma)
    k4_r = martwi(i + h * k3_i, gamma)

    s_nowe = s + (h / 6) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s)
    e_nowe = e + (h / 6) * (k1_e + 2 * k2_e + 2 * k3_e + k4_e)
    i_nowe = i + (h / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i)
    r_nowe = r + (h / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)

    return s_nowe, e_nowe, i_nowe, r_nowe

def model_seir(beta, sigma, gamma, tytuł):
    t_max = 50
    h = 0.1
    kroki = int(t_max / h) + 1

    czas = np.linspace(0, t_max, kroki)
    podatni = np.zeros(kroki)
    narażeni = np.zeros(kroki)
    zainfekowani = np.zeros(kroki)
    martwi = np.zeros(kroki)

    podatni[0] = s_0
    narażeni[0] = e_0
    zainfekowani[0] = i_0
    martwi[0] = r_0

    for i in range(1, kroki):
        podatni[i], narażeni[i], zainfekowani[i], martwi[i] = metoda_rungego_kutty(
            podatni[i - 1], narażeni[i - 1], zainfekowani[i - 1], martwi[i - 1], beta, sigma, gamma, h
        )

    plt.plot(czas, podatni, label="Podatni", color='royalblue')
    plt.plot(czas, narażeni, label="Narażeni", color='lightpink')
    plt.plot(czas, zainfekowani, label="Zainfekowani", color='indianred')
    plt.plot(czas, martwi, label="Martwi", color='mediumvioletred')
    plt.xlabel("Czas")
    plt.ylabel("Populacja")
    plt.legend()
    plt.title(tytuł)
    plt.show()

R0_niskie = 0.5
beta_niskie = R0_niskie * gamma
model_seir(beta_niskie, sigma, gamma, f'Model SEIR dla R0 = {R0_niskie}')

R0_wysokie = 50
beta_wysokie = R0_wysokie * gamma
model_seir(beta_wysokie, sigma, gamma, f'Model SEIR dla R0 = {R0_wysokie}')


'''
R0 ustawione na 0.5 oznacza, że przeciętnie każda zainfekowana osoba zaraża tylko pół innych. 
Gdy R0<1, epidemia nie rozwija się i wygasa.


Ustawienie R0 na 50 wskazuje na znacznie szybsze rozprzestrzenianie się epidemii. 
Choroba silnie wpływa na populację, co objawia się dużą liczbą zainfekowanych oraz 
osób narażonych na działanie wirusa. Wysoka wartość R0 sugeruje zwiększoną zdolność do przenoszenia choroby, 
co prowadzi do bardziej rozległej epidemii.
'''
