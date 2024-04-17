import numpy as np
import matplotlib.pyplot as plt


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


def orden_numpy(n):
    float_str = str(np.format_float_scientific(n, precision=1, exp_digits=2))
    
    return int(float_str.split("e")[1])

def redondear_valor_y_error( valor, error):
    orden_error = orden_numpy(error)
    error_redondeado_hacia_arriba = round(error + 0.5*pow(10,orden_error), abs( orden_error ) )
    
    valor_redondeado = round(valor, -orden_error)
    resultado = [valor_redondeado, error_redondeado_hacia_arriba]
    return resultado


def g(u):
    return ( 2* (np.sinh(u * round(x1, 1)) / u) ) - round(L, 1)


def newton_raphson(u, x1, L):
    N1 = u*x1
    N2 = np.sinh(N1)
    N3 = u*u
    N4 = L*N3
    N5 = N2*u
    N6 = 2*N5
    N7 = np.cosh(N1)
    N8 = N7*N1
    N9 = N8-N2
    N10 = 2*N9
    N11 = N6-N4
    N12 = N11/N10
    N13 = u - N12

    return N13

def buscar_raiz_newton_raphson(u):
    iteraciones = 0
    while iteraciones < 100 and g(u) != 0:
        u_previo = u
        u = newton_raphson(u, x1, L)
        print(u)
        iteraciones = iteraciones + 1
    print("iteraciones:", iteraciones)
    print("Resultado de N-R: ", u)
    return u, u_previo


def error_u_relativo(u, x1, L, ix, iL):
    N1 = u*x1
    N2 = np.sinh(N1)
    N3 = u*u
    N4 = L*N3
    N5 = N2*u
    N6 = 2*N5
    N7 = np.cosh(N1)
    N8 = N7*N1
    N9 = N8-N2
    N10 = 2*N9
    N11 = N6-N4
    N12 = N11/N10
    N13 = u - N12

    p1 = np.cosh(u*x1)*pow(u,2)*x1
    p2 = np.sinh(u*x1)*pow(u,2)*pow(x1,2)
    p3 = np.sinh(u*x1)*u
    p4 = L*pow(u,2)
    p5 = np.cosh(u*x1)*u*x1
    
    coef_x = abs(2*p1/N11-p2/N9)*ix
    coef_L = abs(p4/N11)*iL

    term_1 = (2*abs(p1)+2*abs(p3)*3-2*abs(p4))/abs(N11)
    term_2 = (abs(p2)+abs(p5)*2-abs(np.sinh(u*x1)))/abs(N9)

    em = pow(10,-15)*0.5

    i12 = coef_x*ix+coef_L*iL+(term_1+term_2+4)*em
    error_final = abs(N12*i12)/abs(u-N12)
    return error_final


def forma_analitica():
    semilla = 0.01
    u, u_previo = buscar_raiz_newton_raphson(semilla)
    
    iu_rel = error_u_relativo(u_previo, x1, L, error_x1/x1, error_L/L)
    #busco el error absoluto
    eu = iu_rel*u
    print(u, eu)

    u, eu = redondear_valor_y_error(u, eu)
    print(u, eu)

# AHORA CON PERTURBACIONES


# error de u = error experimental inherente + error de representación
def error_experimental_u_inherente_absoluto(u, e_x1,e_L):
    #e_x1,e_L son las perturbaciones
    coef_x1=(newton_raphson(u, x1 + e_x1, L)-newton_raphson(u, x1 - e_x1, L) )/(2*e_x1)
    coef_L=(newton_raphson(u,x1, L+e_L)-newton_raphson(u,x1, L-e_L) )/(2*e_L)

    coef_x1 = abs(coef_x1)
    coef_L = abs(coef_L)
    #bien redondeado, por eso 0.05
    return coef_x1*0.05 + coef_L*0.05


def perturbacion_te_relativa(u,x1,L):
    ys=newton_raphson( np.float32(u),np.float32(x1),np.float32(L))
    yd=newton_raphson( np.float64(u),np.float64(x1),np.float64(L))
    return abs(yd-ys)/(abs(yd)*0.5*pow(10,-15))


def forma_perturbaciones_experimentales():
    semilla = 0.01
    u, u_previo = buscar_raiz_newton_raphson(semilla)
    
    for i in range(0,3):
        error_exp = error_experimental_u_inherente_absoluto(u_previo, pow(10,-i), pow(10,-i))
        print("variacion de ",pow(10,-i), " error:", error_exp)
    
    #con los datos anteriores determino una cota para el coeficiente del error inherente propagado
    error_exp_inherente = 0.0003
    error_exp_inherente

    # debo manda el u_previo porque estoy propagando la funcion h(u) newton-raphson la cual me da u=h(u_previo)
    te_u = perturbacion_te_relativa(u_previo,x1,L)
    te_u

    error_de_operacion_rel = te_u*0.5*pow(10,-15)
    error_de_operacion_rel

    error_de_operacion_abs = error_de_operacion_rel*newton_raphson(0.045878023115498545,x1,L)
    error_de_operacion_abs

    # Relacion entre el error de propagacion inherente y el error absoluto de redondeo
    error_exp_inherente/error_de_operacion_abs

    #podemos despreciar el error de representación o error de operación.
    u, eu = redondear_valor_y_error(u, error_exp_inherente)
    print(u, eu)


print("DESARROLLO ANALITICO")
forma_analitica()

print("DESARROLLO CON PERTURBACIONES EXPERIMENTALES")
forma_perturbaciones_experimentales()