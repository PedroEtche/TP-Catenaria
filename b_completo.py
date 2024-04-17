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

def mostrar_intervalo_en_g(intervalo, x1, L):
    print(fa.g(intervalo, x1, L))
    print(fa.g(intervalo, x1, L) < 0)

def obtener_cambio_de_signo_en_g(intervalo, x1, L):
    a = 0
    b = 0
    for i, imagen in enumerate(fa.g(intervalo, x1, L)):
        if imagen < 0:
            a = intervalo[i]
            b = intervalo[i+1]
    return (a, b)

def mostrar_valores_de_cambio_de_signo_en_g(a, b):
    print('La funcion en a vale = ', fa.g(a, x1, L))
    print('La funcion en b vale = ', fa.g(b, x1, L))


def graficar_g_en_intervalo(intervalo, x1, L):
    plt.plot(intervalo, fa.g(intervalo, x1, L))
    plt.show()


def verificacion_numerica_N_R(y):
    if np.all( 0.01 < y) and np.all( y < 1):
        print("Imagen dentro del intervalo, se puede usar Newton-Raphson")
    else:
        print("No se puede usar Newton-Raphson")
    

def verificacion_visual_N_R(x, y):
    plt.plot(x,y)
    plt.show()

print("Comprobamos si hay valores negativos en la función g(x) para saber si se puede usar el método de bisección")
intervalo = np.linspace(-1, 1, 100)
mostrar_intervalo_en_g(intervalo, x1, L)
a, b = obtener_cambio_de_signo_en_g(intervalo, x1, L)
mostrar_valores_de_cambio_de_signo_en_g(a, b)
graficar_g_en_intervalo(intervalo, x1, L)    

print("se puede usar el metodo de biseccion")
a=0.01
b=1
fa.metodo_biseccion_segun_cota(a, b, 0.00001, x1, L)

print("probamos si se puede usar el metodo de n-r con un mejor intervalo debido al metodo de biseccion")
x = np.linspace(0.045, 0.046, 100)
x = np.float64(x)
y = fa.newton_raphson(x,x1,L)
plt.plot(x,y)
plt.show()



if np.all( a < y) and np.all( y < b):
    print("La imagen esta dentro del intervalo")
else:
    print("no se puede usar Newton-Raphson")

np.all(fa.newton_raphson_coef_k(x,x1,L)<40)
coef_k = 40
np.all(fa.newton_raphson_coef_d(x,x1,L)<0.009)
coef_d = 0.009
####
# 40*0.009 = 0.36

print("se cumple que es continua, la imagen pertenece al intervalo, y Kd≤1 por lo que podemos usar N-R")
    
semilla = 0.04501
u, u_previo = fa.buscar_raiz_newton_raphson(semilla, x1, L)


iu_rel = fa.error_u_relativo(u_previo, x1, L, 0, 0)
print("Resultado de error relativo de u: ", iu_rel)


eu = iu_rel*u
u, eu = fa.redondear_valor_y_error(u, eu)
print("u:",u," error u:", eu)

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
print("Error absoluto c2: ", ec2)

c2, ec2 = fa.redondear_valor_y_error(c2, ec2)

# Calculo del valor de y(0)
ordenada = fa.fun_catenaria(0,u,c2)
error_ordenada_rel = fa.error_y(0,u,c2,eu/u,0,ec2/c2)
error_ordenada = error_ordenada_rel*ordenada
ordenada, error_ordenada = fa.redondear_valor_y_error(ordenada, error_ordenada)
print("y(0):", ordenada,"error y(0):", error_ordenada)