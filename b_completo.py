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


intervalo = np.linspace(-1, 1, 100)
print(fa.g(intervalo,x1,L))
print(fa.g(intervalo,x1,L) < 0)
a = 0
b = 0
for i, imagen in enumerate(fa.g(intervalo,x1,L)):
    if imagen < 0:
        a = intervalo[i]
        b = intervalo[i+1]

print('La funcion en a vale = ', fa.g(a,x1,L))
print('La funcion en b vale = ', fa.g(b,x1,L))


plt.plot(intervalo, fa.g(intervalo,x1,L))
plt.show()

fa.metodo_biseccion_segun_cota(0.01, 1, 0.00001, x1, L)

a = 0.044
b = 0.046
x = np.linspace(a, b, 100)
x = np.float64(x)
y = fa.newton_raphson(x,x1,L)
plt.plot(x,y)
plt.show()

if np.all( 0.01 < y) and np.all( y < 1):
    print("imagen dentro del intervalo")
else:
    print("no se puede usar Newton-Raphson")

np.all(fa.newton_raphson_coef_k(x,x1,L) < 42)

np.all(fa.newton_raphson_coef_d(x,x1,L)<0.009)

print("Kd = ", 42*0.009, " ≤ 1")
print("se cumple que es continua, la imagen pertenece al intervalo, y Kd ≤ 1")


intervalo = np.linspace(0.01, 1, 100)
intervalo = np.float64(intervalo)
y = fa.newton_raphson(intervalo, x1, L)

if np.all( 0.01 < y) and np.all( y < 1):
    print("Imagen dentro del intervalo, se puede usar Newton-Raphson")
else:
    print("No se puede usar Newton-Raphson")

semilla = 0.04501
u, u_previo = fa.buscar_raiz_newton_raphson(semilla,x1,L)


iu_rel = fa.error_u_relativo(u_previo, x1, L, 0, 0)
print("Resultado de error relativo de u: ", iu_rel)


eu = iu_rel*u
u, eu = fa.redondear_valor_y_error(u, eu)
print( u, eu)

plt.plot(intervalo,y)
plt.show()
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
