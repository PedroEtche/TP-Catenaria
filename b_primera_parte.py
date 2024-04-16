import numpy as np
import matplotlib.pyplot as plt

# CONSTANTES

P = (107241 + 101906)/2

x0 = -P/3500
x1 = P/3500
y0 = P/2100
y1 = y0
L = P/1300

# FUNCIONES

def g(u):
    return ( 2* (np.sinh(u * round(x1, 1)) / u) ) - round(L, 1)
    
def mostrar_intervalo_en_g(intervalo):
    print(g(intervalo))
    print(g(intervalo) < 0)

def obtener_cambio_de_signo_en_g(intervalo):
    a = 0
    b = 0
    for i, imagen in enumerate(g(intervalo)):
        if imagen < 0:
            a = intervalo[i]
            b = intervalo[i+1]
    return (a, b)

def mostrar_valores_de_cambio_de_signo_en_g(a, b):
    print('La funcion en a vale = ', g(a))
    print('La funcion en b vale = ', g(b))


def graficar_g_en_intervalo(intervalo):
    plt.plot(intervalo, g(intervalo))
    plt.show()

# MAIN  

def main():
    intervalo = np.linspace(-1, 1, 100)
    mostrar_intervalo_en_g(intervalo)
    a, b = obtener_cambio_de_signo_en_g(intervalo)
    mostrar_valores_de_cambio_de_signo_en_g(a, b)
    graficar_g_en_intervalo(intervalo)

main()