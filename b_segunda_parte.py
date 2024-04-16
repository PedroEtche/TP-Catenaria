import numpy as np
import matplotlib.pyplot as plt


# CONSTANTES


P = (107241 + 101906)/2
x0 = -29.9
x1 = -x0
y1 = 49.8
y0 = y1
L = 80.4


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


def verificacion_numerica_N_R(y):
    if np.all( 0.01 < y) and np.all( y < 1):
        print("Imagen dentro del intervalo, se puede usar Newton-Raphson")
    else:
        print("No se puede usar Newton-Raphson")
    

def verificacion_visual_N_R(x, y):
    plt.plot(x,y)
    plt.show()


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
    error_final = abs(N13*i12)/abs(u-N13)
    return error_final


def valor_de_c2(x1, y1, u):
    n1 = u * x1
    n2 = np.cosh(n1)
    n3 = n2/u
    return y1 - n3


def fun_catenaria(x,u,c2):
    n1 = u * x
    n2 = np.cosh(n1)
    n3 = n2/u
    return n3 + c2


def fun_catenaria(x,u,c2):
    n1 = u * x
    n2 = np.cosh(n1)
    n3 = n2/u
    return n3 + c2


def error_c2_relativo(u, x1, y, iu, ix1, iy):
    
    p1 = np.sinh(u*x1)*u*x1
    p2 = np.cosh(u*x1)
    p3 = p2/u
    d1 = u*y-np.cosh(u*x1)
    
    coef_u = abs((-p1+p2)/d1)
    coef_x = abs(p1/d1)
    coef_y = abs(y/(y-p3))
    coef_error_representacion = 1+(p1+2*p2)/d1
    coef_error_representacion = abs(coef_error_representacion)

    em = 0.5*pow(10,-15)
    term_inh = coef_u*iu+coef_x*ix1+coef_y*iy
    
    print("razon inherente/representacion:", term_inh/(coef_error_representacion*em))
    return term_inh + coef_error_representacion*em


def error_y(x,u,c2,iu,ix,ic2):
    p1 = abs(np.sinh(u*x)*u*x)
    p2 = abs(np.cosh(u*x))
    p3 = abs(p2/u)
    d1 = abs(np.cosh(u*x)+u*c2)

    em = 0.5*pow(10,-15)
    coef_u = (p1-p2)/d1
    coef_c2 = c2/(c2+p3)
    coef_x = p1/d1

    term = 1+(p1+2*p2)/d1
    return coef_u*iu + coef_c2*ic2 + coef_x*ix + term*em


def main():
    intervalo = np.linspace(0.01, 1, 100)
    intervalo = np.float64(intervalo)
    y = newton_raphson(intervalo, x1, L)

    verificacion_numerica_N_R(y)
    
    semilla = 0.01
    u, u_previo = buscar_raiz_newton_raphson(semilla)


    iu_rel = error_u_relativo(u_previo, x1, L, 0, 0)
    print("Resultado de error relativo de u: ", iu_rel)


    eu = iu_rel*u
    u, eu = redondear_valor_y_error(u, eu)
    print( u, eu)

    verificacion_visual_N_R(intervalo, y)

    # Calculos para c2

    c2 = valor_de_c2(x1,y1,u)
    print("El valor del c2 es: ", c2)

    # Verifico si el valor es correcto evaluando la función de la catenaria en x0

    y0_mediante_funcion = fun_catenaria(x0,u,c2)
    print("Valor de y0, segun el c2 obtenido: ", y0_mediante_funcion)

    # Estimacion de cota de error de c2
    ic2 = error_c2_relativo(u, x1, y0, eu/u, 0, 0)
    print("Error inherente de c2: ", ic2)

    ec2 = c2*ic2
    print("Error de no se que: ", ec2)

    c2, ec2 = redondear_valor_y_error(c2, ec2)

    # Calculo del valor de y(0)
    ordenada = fun_catenaria(0,u,c2)
    error_ordenada_rel = error_y(0,u,c2,eu/u,0,ec2/c2)
    error_ordenada = error_ordenada_rel*ordenada
    ordenada, error_ordenada = redondear_valor_y_error(ordenada, error_ordenada)
    print(ordenada, error_ordenada)

main()