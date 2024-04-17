import numpy as np
import matplotlib.pyplot as plt
import funciones_auxiliares as fa

# CONSTANTES


P = (107241 + 101906)/2
x0 = -29.9
x1 = -x0
y1 = 49.8
y0 = y1
L = 80.4


# FUNCIONES

def mostrar_intervalo_en_g(intervalo):
    print(fa.g(intervalo))
    print(fa.g(intervalo) < 0)

def obtener_cambio_de_signo_en_g(intervalo):
    a = 0
    b = 0
    for i, imagen in enumerate(fa.g(intervalo)):
        if imagen < 0:
            a = intervalo[i]
            b = intervalo[i+1]
    return (a, b)

def mostrar_valores_de_cambio_de_signo_en_g(a, b):
    print('La funcion en a vale = ', fa.g(a))
    print('La funcion en b vale = ', fa.g(b))


def graficar_g_en_intervalo(intervalo):
    plt.plot(intervalo, fa.g(intervalo))
    plt.show()


def verificacion_numerica_N_R(y):
    if np.all( 0.01 < y) and np.all( y < 1):
        print("Imagen dentro del intervalo, se puede usar Newton-Raphson")
    else:
        print("No se puede usar Newton-Raphson")
    

def verificacion_visual_N_R(x, y):
    plt.plot(x,y)
    plt.show()


def main():

    intervalo = np.linspace(-1, 1, 100)
    mostrar_intervalo_en_g(intervalo)
    a, b = obtener_cambio_de_signo_en_g(intervalo)
    mostrar_valores_de_cambio_de_signo_en_g(a, b)
    graficar_g_en_intervalo(intervalo)    
    
    
    
    intervalo = np.linspace(0.01, 1, 100)
    intervalo = np.float64(intervalo)
    y = fa.newton_raphson(intervalo, x1, L)

    verificacion_numerica_N_R(y)
    
    semilla = 0.01
    u, u_previo = fa.buscar_raiz_newton_raphson(semilla)


    iu_rel = fa.error_u_relativo(u_previo, x1, L, 0, 0)
    print("Resultado de error relativo de u: ", iu_rel)


    eu = iu_rel*u
    u, eu = fa.redondear_valor_y_error(u, eu)
    print( u, eu)

    verificacion_visual_N_R(intervalo, y)

    # Calculos para c2

    c2 = fa.valor_de_c2(x1,y1,u)
    print("El valor del c2 es: ", c2)

    # Verifico si el valor es correcto evaluando la función de la catenaria en x0

    y0_mediante_funcion = fa.fun_catenaria(x0,u,c2)
    print("Valor de y0, segun el c2 obtenido: ", y0_mediante_funcion)

    # Estimacion de cota de error de c2
    ic2 = fa.error_c2_relativo(u, x1, y0, eu/u, 0, 0)
    print("Error inherente de c2: ", ic2)

    ec2 = c2*ic2
    print("Error de no se que: ", ec2)

    c2, ec2 = fa.redondear_valor_y_error(c2, ec2)

    # Calculo del valor de y(0)
    ordenada = fa.fun_catenaria(0,u,c2)
    error_ordenada_rel = fa.error_y(0,u,c2,eu/u,0,ec2/c2)
    error_ordenada = error_ordenada_rel*ordenada
    ordenada, error_ordenada = fa.redondear_valor_y_error(ordenada, error_ordenada)
    print(ordenada, error_ordenada)

main()