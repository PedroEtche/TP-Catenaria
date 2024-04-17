import numpy as np
import funciones_auxiliares as fa


# CONSTANTES

x0 = -29.9
x1 = -x0
y1 = 49.8
y0 = y1
L = 80.4

error_x1 = 0.05
error_y1 = 0.05
error_L = 0.05

# FUNCIONES

semilla = 0.04501
u, u_previo = fa.buscar_raiz_newton_raphson(semilla,x1,L)
iu_rel = fa.error_u_relativo(u_previo, x1, L, error_x1/x1, error_L/L)
eu = iu_rel*u
u, eu = fa.redondear_valor_y_error(u, eu)
print("u:",u," error u:", eu)

c2 = fa.fun_c2(x1, y1, u)

error_c2_rel = fa.error_c2_relativo(u, x1, y0, eu/u, 0.05/x1, 0.05/y1)


error_c2 = error_c2_rel*c2
error_c2

c2, error_c2 = fa.redondear_valor_y_error(c2, error_c2)
print("c2 y su error son: ", c2, error_c2)

# Calculos y(0)

ic2 = fa.error_c2_relativo(u, x1, y0, eu/u, 0, 0)
ec2 = c2*ic2
c2, ec2 = fa.redondear_valor_y_error(c2, ec2)

fa.fun_catenaria(0,u,c2)

ordenada = fa.fun_catenaria(0,u,c2)

error_ordenada_rel = fa.error_y(0,u,c2,0,eu/c2,ec2/c2)

error_ordenada = error_ordenada_rel*ordenada

ordenada, error_ordenada = fa.redondear_valor_y_error(ordenada, error_ordenada)
print("La ordenada y su error son: ", ordenada, error_ordenada)