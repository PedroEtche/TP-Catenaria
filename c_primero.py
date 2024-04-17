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

# error de u = error experimental inherente + error de representación
def error_experimental_u_inherente_absoluto(u, e_x1,e_L):
    #e_x1,e_L son las perturbaciones
    coef_x1=(fa.newton_raphson(u, x1 + e_x1, L)-fa.newton_raphson(u, x1 - e_x1, L) )/(2*e_x1)
    coef_L=(fa.newton_raphson(u,x1, L+e_L)-fa.newton_raphson(u,x1, L-e_L) )/(2*e_L)

    coef_x1 = abs(coef_x1)
    coef_L = abs(coef_L)
    #bien redondeado, por eso 0.05
    return coef_x1*0.05 + coef_L*0.05


def perturbacion_te_relativa(u,x1,L):
    ys=fa.newton_raphson( np.float32(u),np.float32(x1),np.float32(L))
    yd=fa.newton_raphson( np.float64(u),np.float64(x1),np.float64(L))
    return abs(yd-ys)/(abs(yd)*0.5*pow(10,-15))


print("DESARROLLO ANALITICO")
semilla = 0.04501
u, u_previo = fa.buscar_raiz_newton_raphson(semilla,x1,L)
    
iu_rel = fa.error_u_relativo(u_previo, x1, L, error_x1/x1, error_L/L)
#busco el error absoluto
eu = iu_rel*u

u, eu = fa.redondear_valor_y_error(u, eu)
print("u:",u," error u:", eu)

print("DESARROLLO CON PERTURBACIONES EXPERIMENTALES")
semilla = 0.01
u, u_previo = fa.buscar_raiz_newton_raphson(semilla,x1,L)

print("Perturbamos los valores de a 1,0.1,0.01")
for i in range(0,3):
    error_exp = error_experimental_u_inherente_absoluto(u_previo, pow(10,-i), pow(10,-i))
    print("variacion de ",pow(10,-i), " error:", error_exp)
    
#con los datos anteriores determino una cota para el coeficiente del error inherente propagado

error_exp_inherente = 0.0003
error_exp_inherente
print("una cota es:0.0003")
# debo manda el u_previo porque estoy propagando la funcion h(u) newton-raphson la cual me da u=h(u_previo)
te_u = perturbacion_te_relativa(u_previo,x1,L)

error_de_operacion_rel = te_u*0.5*pow(10,-15)
    
error_de_operacion_abs = error_de_operacion_rel*fa.newton_raphson(0.045878023115498545,x1,L)

# Relacion entre el error de propagacion inherente y el error absoluto de redondeo
error_exp_inherente/error_de_operacion_abs

#podemos despreciar el error de representación o error de operación.
u, eu = fa.redondear_valor_y_error(u, error_exp_inherente)
print(u, eu)