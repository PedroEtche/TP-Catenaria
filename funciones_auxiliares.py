import numpy as np

error_x1 = 0.05
error_y1 = 0.05
error_L = 0.05

def g(u,x1,L):
    y = ( 2* (np.sinh(u*x1)/u) ) - L
    return y

def orden_numpy(n):
    float_str = str(np.format_float_scientific(n, precision=1, exp_digits=2))
    
    return int(float_str.split("e")[1])

def redondear_valor_y_error( valor, error):
    orden_error = orden_numpy(error)
    error_redondeado_hacia_arriba = round(error + 0.5*pow(10,orden_error), abs( orden_error ) )
    
    valor_redondeado = round(valor, -orden_error)
    resultado = [valor_redondeado, error_redondeado_hacia_arriba]
    return resultado

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

def buscar_raiz_newton_raphson(u, x1, L):
    iteraciones = 0
    while iteraciones < 100 and g(u,x1,L) != 0:
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

def fun_c2(x1, y1, u):
    n1 = u * x1
    n2 = np.cosh(n1)
    n3 = n2/u
    return y1 - n3

def g_derivada(u,x1,L):
    return (2*np.cosh(u*x1)*(u*x1)-2*np.sinh(u*x1))/pow(u,2)

def g_segunda_derivada(u,x1,L):
    j = np.cosh(u*x1)*(u*x1)-np.sinh(u*x1)
    j_der = np.sinh(u*x1)*(u*x1)

    num = j_der*u-2*j
    return (2*num)/pow(u,3)

def newton_raphson_coef_k(x,x1,L):
    return abs(g_segunda_derivada(x,x1,L))/abs(g_derivada(x,x1,L))

def newton_raphson_coef_d(x,x1,L):
    return abs(g(x,x1,L))/abs(g_derivada(x,x1,L))

def metodo_biseccion_segun_cota(a, b, cota, x1, L):

    iteraciones = 0
    while abs(b - a)/2 > cota:
        iteraciones = iteraciones + 1
        c = (a + b) / 2
        if g(a, x1, L) * g(c, x1, L) < 0:
            b = c
        else:
            a = c
        print("a:",a,"g(a):",g(a, x1, L),"  b:",b,"g(b)",g(b, x1, L))
    raiz = (a + b) / 2
    print(f"La raíz de la función es: {raiz:.5f}")
    print("Las iteraciones fueron:", iteraciones)
    return raiz

def funcion_cuadratica(x,a,b,c):
    return a+b*x+c*x*x

def aproximacion_cuadratica(x, y):
    sum = 0
    m = len(x)
    A = [ [1,             x.mean(),      (x*x).mean()],
               [x.mean(),     (x*x).mean(),   (x*x*x).mean()],
               [(x*x).mean(), (x*x*x).mean(), (x*x*x*x).mean()]
             ]
    
    B = [y.mean(), (y*x).mean(), (y*x*x).mean()]
    A = np.array(A)
    B = np.array(B)
    
    X = np.linalg.solve(A, B)
    sol = A*X
    print(sol)
    print(B)
    return X

def error_cuadratico_aproximacion_cuadratica(x,y,a,b,c):
    diferencia = y - funcion_cuadratica(x,a,b,c)
    return np.sqrt( (diferencia*diferencia).mean())

def error_cuadratico_funcion_catenaria(x,y,u,c2):
    diferencia = y - fun_catenaria(x,u,c2)
    return np.sqrt( (diferencia*diferencia).mean())