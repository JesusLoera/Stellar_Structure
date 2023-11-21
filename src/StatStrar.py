# Codigo de starstar.f en Python

# Importamos librerias
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Empieza el programa principal
print("""
¡Bienvenido a StatStar en Python!
      
Ingrese los datos de una estrella para calcular su estructura estelar.    
""")

# Caracteristicas de la estrella
Msolar = float(input("Ingrese la masa (en masas solares): "))
Lsolar = float(input("Ingrese la luminosidad (en luminosidades solares): "))
Te = float(input("Ingrese la temperatura efectiva (en Kelvin): "))
X = float(input("Ingrese la fraccion de masa de hidrogeno (X): "))
Z = float(input("Ingrese la fraccion de metales Z: "))

# Calculamos la fraccion de Helio
Y = 1 - X - Z
while(Y<0):
    print("Las fracciones de H, He y metales deben sumar 1...")
    X = float(input("Ingrese la fraccion de masa hidrogeno X: "))
    Z = float(input("Ingrese la fraccion de metales Z: "))
    Y = 1 - X - Z

# Definicion de constantes
Rsun, Msun, Lsun = 6.9599e10, 1.989e33, 3.826e33 
sigma, c, a =  5.67051e-5, 2.99792458e10, 7.56591e-15
G, k_B, m_H = 6.67259e-8, 1.380658e-16, 1.673534e-24
pi, gamma = 3.141592654e0, 5.0e0/3
tog_bf0, g_ff = 0.01, 1.0e0

# GLOSARIO DE VARIABLES

# deltar = paso de integracion para el radio
# idrflg = flag que define el tamaño del paso de integracion
#        = 0 (tamaño de paso inical para la superficie de Rs/1000.)
#        = 1 (tamaño de paso estandar de Rs/100.)
#        = 2 (tamaño de paso para el nucleo de Rs/5000.)
#  
# Nstart = numero de pasos para los que deben utilizarse ecuaciones de partida
#          (Se supone que la zona mas externa es radiativa)
# Nstop = Numero maximo de zonas permitidas en la estrella.
# Igoof = flag con la condicion final del modelo
#       = -1 (numero de zonas excedido; tambien el valor inicial)
#       =  0 (buen modelo)
#       =  1 (la densidad del nucleo es extrema)
#       =  2 (la luminosidad del nucleo es extrema)
#       =  3 (la temperatura extrapolada del nucleo es muy baja)
#       =  4 (la masa se volvio negativa antes de alcanzar el nucleo)
#       =  5 (la luminosidad se volcio negativa antes de alcanzar el nucleo)
# T0,P0 = temperatura y presion en la superficie (se asume T0 = P0 = 0)


# Parametros
nsh = 999
Nstart, Nstop, Igoof, ierr = 10, nsh, -1, 0
P0, T0, dlPlim, debug = 0.0, 0.0, 99.9, 0 

# Conversion de unidades y calculo del radio solar
# Ms,Ls,Rs = Masa, luminosidad y radio de la estrella (cgs)
Ms = Msolar*Msun
Ls = Lsolar*Lsun
Rs = np.sqrt(Ls/(4.e0*pi*sigma))/(Te**2)
Rsolar = Rs/Rsun

# Definimos funciones que impriman los mensajes de error
def formato200():
    formato = """La variacion de la masa se ha vuelto mas grande de 0.001*Ms dejando el ciclo de aproximacion antes de alcanzar Nstart"""
    print(formato)

def formato300():
    formato = """Ocurrio un problema en la integracion de Runge Kutta"""
    print(formato)

def formato400(r, rho, M_r, kappa, T, epsilon, P, L_r):
    formato = f"""Los valores de la zona anterior son: :
                  r/Rs    = {r/Rs}
                  rho     = {rho} g/cm**3 
                  M_r/Ms  = {M_r/Ms} 
                  kappa   = {kappa} cm**2/g
                  T       = {T} K
                  epsilon = {epsilon} ergs/g/s
                  P       = {P} dinas/cm**2',/,
                  L_r/Ls  = {L_r/Ls}"""
    print(formato)

def formato1000():
    formato = """Un Modelo de Secuencia Principal Homogeneo"""
    print(formato)

def formato2000(Msolar, Mcrat, Rsolar, Rcrat, Lsolar, Lcrat, Te,
    rhocor, X, Tcore, Y, Pcore, Z, epscor, dlpdlt):
    formato = f"""Condiciones de la superficie:          Condiciones centrales:
    Mtot = {Msolar} Msolar            Mc/Mtot     = {Mcrat}
    Rtot = {Rsolar} Rsolar            Rc/Rtot     = {Rcrat}
    Ltot = {Lsolar} Lsolar            Lc/Ltot     = {Lcrat}
    Teff = {Te} K                     Density     = {rhocor} g/cm**3
    X    = {X}                        Temperatura = {Tcore} K
    Y    = {Y}                        Presion     = {Pcore} dinas/cm**2
    Z    = {Z}                        epsilon     = {epscor} ergs/s/g
                                      dlnP/dlnT   = {dlpdlt}"""
    print(formato)
 
def formato2500(Ms):
    formato = f"""Notas:
    (1) La masa aparece como Qm = 1,0 - M_r/Mtot, donde Mtot = {Ms} g
    (2) Las zonas convectivas se indican con c, las zonas radiativas con r
    (3) dlnP/dlnT puede limitarse a +99,9 o -99,9; si es asi esta etiquetado por *)"""
    print(formato)

def formato4000():
    formato = """Lamento ser portador de malas noticias, pero... Tu modelo tiene algunos problemas"""
    print(formato) 

def formato5000():
    formato = """Se ha superado el numero de capas permitidas."""
    print(formato)

def formato5200(rho):
    formato = f"""La densidad del nucleo parece un poco fuera de lugar, la densidad deberia aumentar suavemente hacia el centro. La densidad de la ultima zona calculada fue = {rho} gm/cm**3"""
    print(formato) 

def formato5300():
    formato = """Parece que necesitaras un gas de neutrones degenerado y relatividad general para resolver este nucleo. ¿Quien crees que soy, Einstein?"""
    print(formato)

def formato5400(epsilon):
    formato = f"""El nucleo epsilon parece un poco fuera de lugar y deberia variar suavemente cerca del centro. El valor calculado para la ultima zona fue eps = {epsilon} ergs/g/s"""
    print(formato)

def formato5500(T):
    formato = f"""Su temperatura central extrapolada es demasiado baja, un ajuste un poco mas preciso deberia ser suficiente.
    El valor calculado para la ultima zona fue T = {T} K"""
    print(formato)

def formato5600():
    formato = """¡Creaste una estrella con un agujero en el centro!"""
    print(formato)

def formato5700():
    formato = """¡Esta estrella tiene una luminosidad central negativa!"""
    print(formato) 

def formato5800():
    formato = """¡Llegaste al centro antes de que se agotaran la masa y/o la luminosidad!"""
    print(formato)

def formato6000():
    formato = """Parece que te estas acercando, sin embargo, todavia hay algunos errores menores."""
    print(formato)

def formato7000():
    formato = """¡FELICIDADES, CREO QUE LO ENCONTRASTE! Sin embargo, asegurese de observar tu modelo con atencion."""
    print(formato)

def formato9000():
    formato = """***** La integracion ha sido completada *****
    El modelo se ha almacenado en starmodl.dat.'"""
    print(formato)  


# INICIALIZAMOS VARIBALES

# Fraccion de masa del ciclo CNO definida a 0.5*Z
XCNO = Z/2.0e0

#  Calcula el peso molecular medio mu asumiendo una ionizacion completa
#  (vea Eq. 10.21).
mu = 1.0e0/(2.0*X + 0.75*Y + 0.5*Z)

#  Calcula el delimitador entre conveccion adiabatica y radiacion
#  (vea Eq. 10.87).
gamrat = gamma/(gamma - 1.0e0)

# Crearemos arreglos de tamaño (999)
N = nsh

# Creamos los arreglos que vamos a usar para almacenar
# las soluciones de la estructura estelar
r, P = np.zeros(N, float), np.zeros(N, float)
M_r, L_r = np.zeros(N, float), np.zeros(N, float)
T, rho = np.zeros(N, float), np.zeros(N, float)
kappa, epsilon = np.zeros(N, float), np.zeros(N, float)
tog_bf, dlPdlT = np.zeros(N, float), np.zeros(N, float)

# GLOSARIO DE VARIABLES 2

# Rsun = radio del Sol
# Msun = masa del Sol
# Lsun = luminosidad del Sol
# sigma = constante de Stefan-Boltzmann
# c = velocidad de la luz en el vacio
# a = 4*sigma/c (constante de presion de radiacion)
# G = constante de gravitacion universal
# K_B = constante de Boltzmann
# m_H = masa del atomo de hidrogeno
# pi = 3.141592654
# gamma = 5/3 (gamma adiabatica para un gas monoatomico)
# gamrat = gamma/(gamma-1)
# kPad = P/T**(gamma/(gamma-1)) (constante de un gas adiabatico)
# tog_bf = bound-free opacity constant (ratio of guillotine to gaunt factors)
# g_ff = free-free opacity gaunt factor (assumed to be unity)

#inicializar variables utilizadas para 4 derivadas de ecuaciones de estructura
f_im1=np.zeros(4,float)
dfdr=np.zeros(4,float)
f_i=np.zeros(4,float)

# DEFINIMOS EOS
# Calcula los valores de la densidad, la opacidad, el 
# guillotine-to-gaunt factor ratio y la razon de generacion 
# de la energia en un radio r. 

def EOS(P, T, izone):

    def formato100(P, T, izone):
        formato = f"""Algo anda un poco mal aqui.
        Me estas pidiendo que me ocupe de una temperatura negativa.
        o una presion negativa. ¡Lo siento pero eso no esta en mi
        contrato! Tendras que volver a intentarlo con otra condiciones
        iniciales.
        Por si sirve de ayuda, detecte el problema en la zona {izone}
        con las siguientes condiciones:
        T = {T} K,
        P = {P} dinas/cm**2 """
        print(formato)

    def formato200(P, T, izone, Prad, Pgas, rho):
        formato = f"""'Lo siento, pero se detecto una densidad negativa,
        Mi rutina de ecuacion de estado esta un poco desconcertada por este
        nuevo sistema fisico que has creado. La presion de radiacion
        Probablemente sea demasiado grande, lo que implica que la estrella es
        inestable. Pruebe algo un poco menos radical la proxima vez.
        Por si sirve de ayuda, detecte el problema en la zona {izone}
        con las siguientes condiciones:
        T       = {T} K ,
        P_total = {P} dinas/cm**2 ,
        P_rad   = {Prad} dinas/cm**2 ,
        P_gas   = {Pgas} dinas/cm**2 ,
        rho     = {rho} g/cm**3 """
        print(formato)


    oneo3=0.333333333e0
    twoo3=0.666666667e0

    if ((T < 0.0e0) or (P < 0.0e0)):
        formato100(P, T, izone)
        return (0.0, 0.0, 0.0, 0.0, 1)

    Prad = a*T**4/3.0e0
    Pgas = P - Prad
    rho = (mu*m_H/k_B)*(Pgas/T)

    if (rho < 0.0e0):
        formato200(P, T, izone, Prad, Pgas, rho)
        return (rho, 0.0 , 0.0 ,0.0, 1)
    
    # Calcular la opacidad, incluida la relacion del factor
    # guillotine-to-gaunt; Novotny (1973), pag. 469. k_bf, k_ff y k_e son las
    # opacidades libres-ligadas, libres-libres y de dispersion de electrones,
    # dadas por las ecuaciones. (9.19), (9.20) y (9.21), respectivamente.

    tog_bf = 2.82e0*(rho*(1.0e0 + X))**0.2e0
    k_bf = 4.34e25/tog_bf*Z*(1.0e0 + X)*rho/T**3.5e0
    k_ff = 3.68e22*g_ff*(1.0e0 - Z)*(1.0e0 + X)*rho/T**3.5e0
    k_e = 0.2e0*(1.0e0 + X)
    kappa = k_bf + k_ff + k_e

    # Calcular la generacion de energia por la cadena pp y el ciclo CNO. Estos se calculan usando las Ecs. (10.49) y (10.53), que provienen de Fowler, Caughlan y Zimmerman (1975). El factor de deteccion para la cadena pp se calcula como fpp; vease Clayton (1968), pag. 359ff.

    T6 = T*1.0e-06
    fx = 0.133e0*X*np.sqrt((3.0e0 + X)*rho)/T6**1.5e0
    fpp = 1.0e0 + fx*X
    psipp = 1.0e0 + 1.412e8*(1.0e0/X - 1.0e0)*np.exp(-49.98*T6**(-oneo3))
    Cpp = 1.0e0 + 0.0123e0*T6**oneo3 + 0.0109e0*T6**twoo3 + 0.000938e0*T6
    epspp = 2.38e6*rho*X*X*fpp*psipp*Cpp*T6**(-twoo3)*np.exp(-33.80e0*T6**(-oneo3))
    CCNO = 1.0e0 + 0.0027e0*T6**oneo3 - 0.00778e0*T6**twoo3 - 0.000149e0*T6
    epsCNO = 8.67e27*rho*X*XCNO*CCNO*T6**(-twoo3)*np.exp(-152.28e0*T6**(-oneo3))
    epsilon = epspp + epsCNO

    return (rho, kappa, epsilon, tog_bf, 0)


# Definimos STARTMDL
def STARTMDL(r_i, M_ri, L_ri, tog_bf, irc):

        r = r_i + deltar
        M_rip1 = M_ri
        L_rip1 = L_ri

        # Esta es la aproximacion radiativa (desprecie la presion de radiacion
        # y la opacidad de dispersion de electrones) ver Prialnik Eq. 5.1, 5.3
        # y Secc. 3.7 o C&O Ecs. (H.1), (H.2), (9.19) y (9.20).

        if (irc == 0):
            T_ip1 = G*M_rip1*mu*m_H/(4.25e0*k_B)*(1.0e0/r - 1.0e0/Rs)
            A_bf = 4.34e25*Z*(1.0e0 + X)/tog_bf
            A_ff = 3.68e22*g_ff*(1.0e0 - Z)*(1.0e0 + X)
            Afac = A_bf + A_ff
            P_ip1 = np.sqrt((1.0e0/4.25e0)*(16.0e0/3.0e0*np.pi*a*c)*(G*M_rip1/L_rip1)*(k_B/(Afac*mu*m_H)))*T_ip1**4.25e0

        # Esta es la aproximacion convectiva
        # ver Prialnik Sec 6.5, 6.6 o C&O Eqns. (H.3) y (10.75).

        else:
            T_ip1 = G*M_rip1*mu*m_H/k_B*(1.0e0/r - 1.0e0/Rs)/gamrat
            P_ip1 = kPad*T_ip1**gamrat

        return r,P_ip1, M_rip1, L_rip1, T_ip1

#  Las 4 funciones calculan los gradientes de presion, masa, luminosidad y temperatura en un radio r.

def dPdr(r, M_r, rho):
    return -G*rho*M_r/r**2
    
def dMdr(r, rho):
    return (4.0e0*np.pi*rho*r**2)
    
def dLdr(r, rho, epsilon):
    return (4.0e0*np.pi*rho*epsilon*r**2)

def dTdr(r, M_r, L_r, T, rho, kappa, irc):
    if (irc == 0):
        return (-(3.0e0/(16.0e0*np.pi*a*c))*kappa*rho/T**3*L_r/r**2)
    # Este es el gradiente de temperatura convectiva adiabatica 
    # (Prialnik Eq. 6.29 or C&O Eq. 10.81).
    else:
        return (-1.0e0/gamrat*G*M_r/r**2*mu*m_H/k_B)

# FUNCIoN FUNDEQ(r, f, irc, izone)
def FUNDEQ(r, f, irc, izone):
    dfdr=np.zeros(4)
    P   = f[0]
    M_r = f[1]
    L_r = f[2]
    T   = f[3]
    rho, kappa, epsilon, tog_bf, ierr = EOS(P, T, izone)
    dfdr[0] = dPdr(r, M_r, rho)
    dfdr[1] = dMdr(r, rho)
    dfdr[2] = dLdr(r, rho, epsilon)
    dfdr[3] = dTdr(r, M_r, L_r, T, rho, kappa, irc)
    return (dfdr,ierr)

# MeTODO DE RUNGE KUTTA 
def RUNGE(f_im1, dfdr, r_im1, deltar, irc, izone):
    f_temp=np.zeros(4)
    f_i=np.zeros(4)
    dr12 = deltar/2.0e0
    dr16 = deltar/6.0e0
    r12  = r_im1 + dr12
    r_i  = r_im1 + deltar
    # Calcular derivadas intermedias a partir de las cuatro estelares
    # fundamentales encontradas en el metodo FUNDEQ.
    for i in range(0,4):
        f_temp[i] = f_im1[i] + dr12*dfdr[i]

    df1, ierr = FUNDEQ(r12, f_temp, irc, izone)
    if (ierr != 0):
        return f_i,ierr

    for i in range(0,4):
        f_temp[i] = f_im1[i] + dr12*df1[i]

    df2, ierr = FUNDEQ(r12, f_temp, irc, izone)

    if (ierr != 0):
        return f_i,ierr

    for i in range(0,4):
        f_temp[i] = f_im1[i] + deltar*df2[i]

    df3, ierr=FUNDEQ(r_i, f_temp, irc, izone)
    if (ierr != 0):
        return f_i,ierr

    # Calcula las variables en la siguiente capa (i + 1).

    for i in range(0,4):
        f_i[i] = f_im1[i] + dr16*(dfdr[i] + 2.0e0*df1[i] + 2.0e0*df2[i] + df3[i])

    return f_i,0


# INTEGRACIoN DE LA ESTRUCTURA ESTELAR

# Empezamos con un muy pequeño paso de integracion debido a que las condiciones en la superficie varian muy rapido
deltar = -Rs/1000.0e0
idrflg = 0
initsh = 0

# Definimos las condiciones de frontera en 
# la superficie de la estrella
r[initsh], M_r[initsh] = Rs, Ms
L_r[initsh], T[initsh] = Ls, T0
P[initsh], tog_bf[initsh] = P0, tog_bf0

if (P0 <= 0.0) or (T0 <= 0.0):
    rho[initsh]    = 0.0
    kappa[initsh]  = 0.0
    epsilon[initsh] = 0.0
    tog_bf[initsh] = 0.01
else:
    rho[initsh], kappa[initsh], epsilon[initsh], tog_bf[initsh], ierr = EOS(P[initsh], T[initsh], initsh)
    if (ierr != 0):
        print ("Vamos a detenernos")
        istop=0

# Aplicar soluciones superficiales aproximadas para comenzar la integracion,
# suponiendo transporte de radiacion en la zona mas externa (el bucle do 20).
# irc = 0 para radiacion, irc = 1 para conveccion.
# Supongamos valores iniciales arbitrarios para kPad y dlPdlT.
# dlPdlT = dlnP/dlnT (consulte la ecuacion de Prialnik 6.28 o la 
# ecuacion de C&O 10.87)

kPad = 0.3e0
irc = 0
dlPdlT[initsh] = 4.25e0

for i in range(0, Nstart):

    ip1 = i + 1

    r[ip1], P[ip1], M_r[ip1], L_r[ip1], T[ip1]= STARTMDL(r[i], M_r[i], L_r[i], tog_bf[i], irc)

    rho[ip1], kappa[ip1], epsilon[ip1], tog_bf[ip1], ierr = EOS(P[ip1], T[ip1], ip1)


    if (ierr != 0):
        formato400(r[i], rho[i], M_r[i], 
                    kappa[i], T[i], epsilon[i],
                    P[i],
                    L_r[i])      
        break                       

            
    # Determinar si la conveccion funcionara en la siguiente zona
    # calculando dlnP/dlnT numericamente entre las zonas i e i+1 [ip1].
    # Actualice la constante adiabatica del gas si es necesario.

    if (i > initsh):
        dlPdlT[ip1] = np.log(P[ip1]/P[i])/np.log(T[ip1]/T[i])
    else:
        dlPdlT[ip1] = dlPdlT[i]

    if (dlPdlT[ip1] < gamrat):
        irc = 1
    else:
        irc = 0
        kPad = P[ip1]/T[ip1]**gamrat

    # Pruebe para ver si la suposicion masa constante en la superficie sigue siendo valido.

    deltaM = deltar*dMdr(r[ip1], rho[ip1])
    M_r[ip1] = M_r[i] + deltaM
    if (np.abs(deltaM) > (0.001e0*Ms)):
        if (ip1 > 1):
            ip1 = ip1 - 1
            print('La variacion de masa ha llegado a ser mayor que 0.001*Ms')
            print('dejando el bucle de aproximacion antes de que se alcanzara Nstart')
            break

Nsrtp1 = ip1 + 1

if (ierr != 0):    
    # salir si hemos llegado a este punto despues de un error en
    # la inicializacion
    Nstop=Nsrtp1-1
    istop=Nstop

for i in range(Nsrtp1, Nstop):
    im1 = i - 1
    f_im1[0] = P[im1]
    f_im1[1] = M_r[im1]
    f_im1[2] = L_r[im1]
    f_im1[3] = T[im1]
    dfdr[0]  = dPdr(r[im1], M_r[im1], rho[im1])
    dfdr[1]  = dMdr(r[im1], rho[im1])
    dfdr[2]  = dLdr(r[im1], rho[im1], epsilon[im1])
    dfdr[3]  = dTdr(r[im1], M_r[im1], L_r[im1], T[im1], rho[im1], kappa[im1], irc)
    f_i, ierr = RUNGE(f_im1, dfdr, r[im1], deltar, irc, i)

    if (ierr != 0):
        formato300()
        formato400(r[im1], rho[im1], M_r[im1], kappa[im1], T[im1], epsilon[im1], P[im1], L_r[im1])
        break

    # Actualizar los parametros estelares para la siguiente zona, incluyendo
    # la adicion de dr al radio antiguo (notese que dr < 0 ya que la
    # integracion es hacia el interior).

    r[i]   = r[im1] + deltar
    P[i]   = f_i[0]
    M_r[i] = f_i[1]
    L_r[i] = f_i[2]
    T[i]   = f_i[3]

    # Calcule la densidad, opacidad y tasa de generacion de energia 
    # para esta zona.

    rho[i], kappa[i], epsilon[i], tog_bf[i], ierr = EOS(P[i], T[i], i)

    if (ierr != 0):
        formato400(r[im1], rho[im1], M_r[im1], kappa[im1], T[im1], epsilon[im1], P[im1], L_r[im1])
        istop = i
        break

    if (debug == 1): 
        print (i, r[i], M_r[i], L_r[i], T[i], P[i], rho[i], kappa[i], epsilon[i], tog_bf[i])

    # Determinar si la conveccion funcionara en la siguiente zona
    # calculando dlnP/dlnT y comparandolo con gamma/(gamma-1)
    # (consulte la ecuacion de Prialnik 6.28 o la ecuacion de C&O 10.87).
    # Configure la bandera de conveccion apropiadamente.

    dlPdlT[i] = np.log(P[i]/P[im1])/np.log(T[i]/T[im1])

    if (dlPdlT[i] < gamrat):
        irc = 1
    else:
        irc = 0

    # Comprueba si se ha alcanzado el centro.  Si es asi, establece Igoof y
    # estimar las condiciones centrales rhocor, epscor, Pcore, y Tcore.
    # La densidad central se estima que es la densidad media de la
    # bola central restante, la presion central se determina
    # utilizando la expansion de Taylor en el centro (Prialnik - Ejercicio. 5.
    # 1; CO Ec. H.4) y el valor central de la tasa de generacion
    # de generacion de energia se calcula como la luminosidad interior restante
    # luminosidad interior dividida por la masa de la bola central. 
    # Finalmente, la temperatura central se calcula aplicando la ley de los
    # gases ideales (despreciando la presion de radiacion).

    if ((r[i] <= np.abs(deltar)) and ((L_r[i] >= (0.1e0*Ls)) or (M_r[i] >= (0.01e0*Ms)))):
        # Llega al centro antes de que se agote la masa/luminosidad
        Igoof = 6
    elif (L_r[i] <= 0.0e0):
        # Obtenido luminosidad central negativa
        Igoof = 5
        rhocor = M_r[i]/(4.0e0/3.0e0*np.pi*r[i]**3)
        if (M_r[i] != 0.0e0):
            # razon de generacion de energia
            epscor = L_r[i]/M_r[i]
        else:
            epscor = 0.0e0
            Pcore = P[i] + 2.0e0/3.0e0*np.pi*G*rhocor**2*r[i]**2
            Tcore = Pcore*mu*m_H/(rhocor*k_B)
    elif (M_r[i] <= 0.0e0):
        # El modelo tiene un hoyo en el centro (densidad negativa!)
        Igoof  = 4  
        rhocor = 0.0e0
        epscor = 0.0e0
        Pcore  = 0.0e0
        Tcore  = 0.0e0
    elif ((r[i] < (0.02e0*Rs)) and ((M_r[i] < (0.01e0*Ms)) and ((L_r[i] < 0.1e0*Ls)))):
        # si hemos alcanzado <2% del radio de la estrella, 
        # <1% de masa encerrada y <10% de luminosidad entonces.... 
        # rho: expansion de Taylor en el centro
        rhocor = M_r[i]/(4./3.*np.pi*r[i]**3)
        # establecer una masa central maxima razonable
        rhomax = 10.0e0*(rho[i]/rho[im1])*rho[i]
        epscor = L_r[i]/M_r[i]
        # P: expansion de Taylor en el centro
        Pcore  = P[i] + 2.0e0/3.0e0*np.pi*G*rhocor**2*r[i]**2
        # Se asume un gas ideal
        Tcore  = Pcore*mu*m_H/(rhocor*k_B)
        # En general, todos estos deberian producir valores que se 
        # eleven hacia el centro (pero no demasiado altos)
        if ((rhocor < rho[i]) or (rhocor > rhomax)):
            # rho esta fuera de lugar ya sea grande o pequeño
            Igoof = 1
        elif (epscor < epsilon[i]):
            # tasa de generacion de energia un poco baja (baja)
            Igoof = 2
        elif (Tcore < T[i]):
            # Temperatura un poco baja (baja)
            Igoof = 3
        else:
            # Se ha excedido el numero de capas permitidas
            Igoof = 0

    if (Igoof != -1):
        istop = i
        break

    # ¿Es hora de cambiar el tamaño del paso?

    if ((idrflg == 0) and (M_r[i] < (0.99e0*Ms))):
        deltar = (-1.0)*Rs/100.0e0
        idrflg = 1

    if ((idrflg == 1) and (deltar >= (0.5*r[i]))):
        deltar = (-1.0)*Rs/5000.0e0
        idrflg = 2

    istop = i

#  Generate warning messages for the central conditions.

rhocor = M_r[istop]/(4.0e0/3.0e0*np.pi*r[istop]**3)
epscor = L_r[istop]/M_r[istop]
Pcore  = P[istop] + 2.0e0/3.0e0*np.pi*G*rhocor**2*r[istop]**2
Tcore  = Pcore*mu*m_H/(rhocor*k_B)

if  (Igoof != 0):

    if (Igoof == -1):
        formato4000()
        formato5000()

    if (Igoof == 1):
        formato6000()
        formato5200(rho[istop])
        print (rhocor,rhomax)

    if (rhocor > 1e10):
        formato5300()

    if (Igoof == 2):
        formato6000()
        formato5400(epsilon[istop])

    if (Igoof == 3):
        formato6000()
        formato5500(T[istop])

    if (Igoof == 4):
        formato4000()
        formato5600()

    if (Igoof == 5):
        formato4000()
        formato5700()

    if (Igoof == 6):
                formato4000()
                formato5800()
else:
    formato7000()

# Imprimir las condiciones centrales. Si es necesario, establezca limites
# para el radio central, la masa y la luminosidad si es necesario, para
# evitar desbordamientos del campo de formato.

Rcrat = r[istop]/Rs
if (Rcrat < -9.999e0): 
    Rcrat = -9.999e0
        
Mcrat = M_r[istop]/Ms
if (Mcrat < -9.999e0): 
    Mcrat = -9.999e0
        
Lcrat = L_r[istop]/Ls
if (Lcrat < -9.999e0): 
    Lcrat = -9.999e0

f=open('statstar.dat', 'w')

f.write('Un Modelo de Secuencia Principal Homogeneo\n')
f.write(' Condiciones de superficie:     Condiciones centrales:\n')
f.write(' Mtot = {0:13.6E} Msolar        Mc/Mtot     = {1:12.5E}\n'.format(Msolar,Mcrat))
f.write(' Rtot = {0:13.6E} Rsolar        Rc/Rtot     = {1:12.5E}\n'.format(Rsolar,Rcrat))
f.write(' Ltot = {0:13.6E} Lsolar        Lc/Ltot     = {1:12.5E}\n'.format(Lsolar,Lcrat))
f.write(' Teff = {0:13.6E} K             Densidad    = {1:12.5E}\n'.format(Te,rhocor))
f.write(' X    = {0:13.6E}               Temperatura = {1:12.5E}\n'.format(X,Tcore))
f.write(' Y    = {0:13.6E}               Presion:    = {1:12.5E} dinas/cm**2\n'.format(Y,Pcore))
f.write(' Z    = {0:13.6E}               epsilon     = {1:12.5E} ergs/s/g\n'.format(Z,epscor))
f.write('                                dlnP/dlnT   = {0:12.5E}\n'.format(dlPdlT[istop]))


f.write('Notas:\n')
f.write(' (1) La masa aparece como Qm = 1,0 - M_r/Mtot, donde Mtot = {0:13.6}\n'.format(Msun))
f.write(' (2) Las zonas convectivas se indican con c, las zonas radiativas con r\n')
f.write(' (3) dlnP/dlnT puede limitarse a +99,9 o -99,9; si es asi esta\n')
f.write(' etiquetado por *\n')

# Imprime datos desde el centro de la estrella hacia afuera, etiquetando zonas
# convectivas o radiativas con c o r, respectivamente. Si abs(dlnP/dlnT) excede
# 99,9, establezca un indicador de advertencia de impresion (*) y establezca el
# limite de salida en +99,9 o -99,9 segun corresponda para evitar
# desbordamientos de campos de formato.

f.write('   r        Qm       L_r       T        P        rho      kap      eps     dlPdlT\n')

for ic in range(0,istop+1):
    i = istop - ic
    Qm = 1.0e0 - M_r[i]/Ms    # Total mass fraction down to radius

    if (dlPdlT[i] < gamrat):
        rcf = 'c'
    else:
                rcf = 'r'
    if (np.abs(dlPdlT[i]) > dlPlim):
        dlPdlT[i] = np.copysign(dlPlim, dlPdlT[i])
        clim = '*'
    else:
        clim = ' '
    s='{0:7.2E} {1:7.2E} {2:7.2E} {3:7.2E} {4:7.2E} {5:7.2E} {6:7.2E} {7:6.2E}{8:1s}{9:1s} {10:5.1f}\n'.format(r[i], Qm, L_r[i], T[i], P[i], rho[i], kappa[i], epsilon[i], clim, rcf, dlPdlT[i])
    f.write(s)

# Output en pantalla

formato9000()

import pandas as pd
import matplotlib.pyplot as plt

def buscar_cols(file):
    with open(file,"r",encoding="utf-8") as file:
        rows = [row for row in file]
        for count, row in enumerate(rows):
            if (('r' in row)   and
                ('Qm' in row)  and 
                ('L_r' in row) and 
                ('T' in row)   and
                ('dlPdlT' in row)):
                count = count +1
                return count
            
columns = ['r', 'Qm', 'L_r', 'T', 'P', 'rho', 'kap', 'eps', 'zone', 'dlPdlT']

archivo = "statstar.dat"
datos = pd.read_csv(archivo, skiprows=buscar_cols(archivo), sep='\s+', names=columns)

# Grafica Luminosidad vs Radio
plt.plot(datos.r/1000, datos['L_r']) 
plt.xlabel('r [km]')
plt.ylabel(r'$L_r$ $[erg/s]$')
plt.title('Luminosidad')
plt.savefig('luminosidad.png')
plt.clf()

# Grafica Temperatura vs Radio
plt.plot(datos.r/1000, datos['T']) 
plt.xlabel('r [km]')
plt.ylabel(r'$T$ $[K]$')
plt.title('Temperatura')
plt.savefig('temperatura.png')
plt.clf()

# Grafica densidad vs Radio
plt.plot(datos.r/1000, datos['rho']) 
plt.xlabel('r [km]')
plt.ylabel(r'$rho$ $[\frac{g}{cm^3}]$')
plt.title('Densidad')
plt.savefig('densidad.png')
plt.clf()

# Grafica Masa vs Radio
plt.plot(datos.r/1000, datos['Qm']) 
plt.xlabel('r [km]')
plt.ylabel(r'$Q_m$')
plt.title('Masa (Qm)')
plt.savefig('masa.png')
plt.clf()

# Grafica Presión vs Radio
plt.plot(datos.r/1000, datos['P']) 
plt.xlabel('r [km]')
plt.ylabel(r'$P$ $[dinas/cm**2]$')
plt.title('Presión')
plt.savefig('presion.png')
plt.clf()

# Grafica Opacidad vs Radio
plt.plot(datos.r/1000, datos['kap']) 
plt.xlabel('r [km]')
plt.ylabel(r'$\kappa$')
plt.title('Opacidad')
plt.savefig('opacidad.png')
plt.clf()

# Grafica Epsilon vs Radio
plt.plot(datos.r/1000, datos['eps']) 
plt.xlabel('r [km]')
plt.ylabel(r'$\epsilon$ ergs/s/g')
plt.title('Epsilon')
plt.savefig('epsilon.png')
plt.clf()

# Grafica dlnP/dlnT vs Radio
plt.plot(datos.r/1000, datos['dlPdlT']) 
plt.xlabel('r [km]')
plt.ylabel(r'dlnP/dlnT')
plt.title('dlnP/dlnT')
plt.savefig('dlnP_dlnT.png')
plt.clf()
