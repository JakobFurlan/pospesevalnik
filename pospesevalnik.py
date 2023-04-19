import numpy as np
import matplotlib.pyplot as plt

# Funkcija za izračun Lorentzove sile na nabit delec
def lorentzova_sila(q, v, B):
    F = q * np.cross(v, B)
    return F

# Vnosi uporabnika
q = float(input("Naboj delca (Coulomb): "))
m = float(input("Masa delca (kg): "))
v0 = np.array([float(input("x komponenta začetna hitrosti (m/s): ")),
              float(input("y komponenta začetne hitrosti (m/s): ")),
              float(input("z komponenta začetne hitrosti (m/s): "))])
r0 = np.array([0.0, 0.0, 0.0])
B = np.array([float(input("x komponenta magnetnega polja (Tesla): ")),
             float(input("y komponenta magnetnega polja (Tesla): ")),
             float(input("z komponenta magnetnega polja (Tesla): "))])
t_max = float(input("celoten čas simulacije (s): "))
dt = float(input("časovni interval v meritev (s): "))
       
# Konstante
c = 2.998e8  # Hitrost svetlobe (m/s)
g = 9.81  #gravitacijski pospešek (m^2/s)

# Začetni podatki
v = v0
r = r0
a = np.zeros(3)
t = 0.0
cas = [t]
pozicija = [r]
hitrost = [v]
pospesek = [a]

# Numerična integracija po Verletovi metodi
while t < t_max:
    # Izračun pospeška pri trenutnem casu
    F = lorentzova_sila(q, v, B)
    a = F / m

    # Posodobitev položaja in hitrosti z Verletovo metodo
    r_naslednji = r + v*dt + 0.5*a*dt**2
    F_naslednji = lorentzova_sila(q, v, B)
    a_naslednji = F_naslednji / m
    v_naslednji = v + 0.5*(a + a_naslednji)*dt
    
    # Preverjanje ali velikost vektorja hitrosti presega svetlobno hitrost
    v_mag = np.linalg.norm(v_naslednji)
    if v_mag >= c:
        # Prilagoditev hitrosti tako, da bo enaka svetlobni
        v_naslednji = v_naslednji / v_mag * c
    r = r_naslednji
    v = v_naslednji
    a = a_naslednji

    # Dodajanje rezultatov k seznamom
    t += dt
    cas.append(t)
    pozicija.append(r)
    hitrost.append(v)
    pospesek.append(a)

# Pretvorba seznamov v numpy nize za lažjo upravljanje s podatki
pozicija = np.array(pozicija)
hitrost = np.array(hitrost)
pospesek = np.array(pospesek)

# Izračun kinetične in potencialne energije
kin_energija = 0.5 * m * (hitrost ** 2).sum(axis=1)
pot_energija = m * g * pozicija[:, 1]

# Izračun skupne energije
skupna_energija = kin_energija + pot_energija

# Izpis rezultatov z grafi
fig, axs = plt.subplots(4, 1, figsize=(10, 10))
axs[0].plot(pozicija[:,0], pozicija[:,1])
axs[0].set_xlabel("x pozicija (m)")
axs[0].set_ylabel("y pozicija (m)")
axs[0].set_title("pozicija")

axs[1].plot(cas, hitrost[:,0], label="x-hitrost")
axs[1].plot(cas, hitrost[:,1], label="y-hitrost")
axs[1].plot(cas, hitrost[:,2], label="z-hitrost")
axs[1].set_xlabel("cas (s)")
axs[1].set_ylabel("hitrost (m/s)")
axs[1].legend()
axs[1].set_title("hitrost")

axs[2].plot(cas, pospesek[:,0], label="x-pospešek")
axs[2].plot(cas, pospesek[:,1], label="y-pospešek")
axs[2].plot(cas, pospesek[:,2], label="z-pospešek")
axs[2].set_xlabel("cas (s)")
axs[2].set_ylabel("pospešek (m/s^2)")
axs[2].legend()
axs[2].set_title("pospešek")

axs[3].plot(cas, kin_energija, label="Kinetična energija")
axs[3].plot(cas, pot_energija, label="Potencialna energija")
axs[3].plot(cas, skupna_energija, label="Skupna energija")
axs[3].set_xlabel("cas (s)")
axs[3].set_ylabel("Energija (J)")
axs[3].legend()
axs[3].set_title("energija")

plt.subplots_adjust(hspace=1)

plt.show()

#priporočene vrednosti
#Naboj delca (Coulomb): 1.6e-19
#Masa delca (kg): 1.67e-27
#x komponenta začetna hitrosti (m/s): 1000000 (1.0e6)
#y komponenta začetne hitrosti (m/s): 0
#z komponenta začetne hitrosti (m/s): 0
#x komponenta magnetnega polja (Tesla): 0
#y komponenta magnetnega polja (Tesla): 0
#z komponenta magnetnega polja (Tesla): 1
#celoten čas simulacije (s): 0.000001 (1.0e-6)
#časovni interval v meritev (s): 0.000000001 (1.0e-9)