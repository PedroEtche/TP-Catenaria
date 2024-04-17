import numpy as np
import matplotlib.pyplot as plt
import funciones_auxiliares as fa

x0 = -12
x1 = -x0
y1 = 66
y0 = y1
L = 120.75

x = np.linspace(-1, 1, 20)
x = np.float64(x)
y = fa.g(x,x1,L)
plt.plot(x,y)
plt.show()

print(y < 0)

a = 0.1
b = 1
x = np.linspace(a, b, 100)
x = np.float64(x)
y = fa.newton_raphson(x,x1,L)
plt.plot(x,y)
plt.show()

if np.all( a < y) and np.all( y < b):
    print("imagen dentro del intervalo")
else:
    print("no se puede usar Newton-Raphson")


np.all(fa.newton_raphson_coef_k(x,x1,L)< 6)

np.all(fa.newton_raphson_coef_d(x,x1,L)<0.01)


# 6*0.01 = 0.06

# se cumple que es continua, la imagen pertenece al intervalo, y Kd≤1


semilla = 0.1
u, u_previo = fa.buscar_raiz_newton_raphson(semilla,x1,L)


iu_rel = fa.error_u_relativo(u_previo, x1, L, 0, 0)

eu = iu_rel*u

u, eu = fa.redondear_valor_y_error(u, eu)
print( u, eu)


c2 = fa.fun_c2(x1,y1,u)

#Verifico si el valor es correcto evaluando la función de la catenaria en x0
y0_mediante_funcion = fa.fun_catenaria(x0,u,c2)
print(y0_mediante_funcion)


ic2 = fa.error_c2_relativo(u, x1, y0, eu/u, 0, 0)
print(ic2)
ec2 = c2*ic2

c2, ec2 = fa.redondear_valor_y_error(c2, ec2)
print(c2, ec2)


ordenada = fa.fun_catenaria(0,u,c2)

error_ordenada_rel = fa.error_y(0,u,c2,eu/u,0,ec2/c2)
error_ordenada = error_ordenada_rel*ordenada
ordenada, error_ordenada = fa.redondear_valor_y_error(ordenada, error_ordenada)
print(ordenada, error_ordenada)  


x_medido = [-12,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,12]
y_medido = [65.212,30.255,18.472,13.56,10.795,9.993,10.441,12.223,17.753,28.369,65.442]
x_medido = np.array(x_medido)
y_medido = np.array(y_medido)
plt.plot(x_medido,y_medido)
plt.show()


a,b,c = fa.aproximacion_cuadratica(x_medido,y_medido)

x = np.linspace(-12,12, 100)
y_cat = fa.fun_catenaria(x,u,c2)
y_cua = fa.funcion_cuadratica(x,a,b,c)

plt.grid(True)
plt.plot(x,y_cat, label='catenaria', color='blue')
plt.plot(x,y_cua, label='cuadratica', color='red')
plt.show()


print(fa.error_cuadratico_aproximacion_cuadratica(x_medido, y_medido, a, b, c))

print(fa.error_cuadratico_funcion_catenaria(x_medido, y_medido, u, c2))