
""" Librerías """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

""" Clase Stellar Structure"""

class StellarStructure():

    def __init__(self, Msolar, Lsolar, Te, X, Z):

 
        """ Constrcutor con los parámetros de la estrella """

        # Características de la estrella
        self.Msolar = Msolar      # Masa (unidades solares)
        self.Lsolar = Lsolar        # Luminosidad (unidades solares)
        self.Te = Te            # Temperatura efectiva (kelvins)
        self.X = X                # Fracción de hidrogeno (adimensional)
        self.Z = Z                # Fracción de metales (adimensional)
        self.Y = 1.0e0 - X - Z    # Fracción de helio (adimensional)

        # deltar = radius integration step
        # idrflg = set size flag
        #        = 0 (initial surface step size of Rs/1000.)
        #        = 1 (standard step size of Rs/100.)
        #        = 2 (core step size of Rs/5000.)
        #  
        # Nstart = number of steps for which starting equations are to be used
        #          (the outermost zone is assumed to be radiative)
        # Nstop = maximum number of allowed zones in the star
        # Igoof = final model condition flag
        #       = -1 (number of zones exceeded; also the initial value)
        #       =  0 (good model)
        #       =  1 (core density was extreme)
        #       =  2 (core luminosity was extreme)
        #       =  3 (extrapolated core temperature is too low)
        #       =  4 (mass became negative before center was reached)
        #       =  5 (luminosity became negative before center was reached)
        # T0, P0 = surface temperature and pressure (T0 = P0 = 0 is assumed)


        """ 
        Definición de constantes en el sistema CGS
        """
        self.Rsun, self.Msun, self.Lsun = 6.9599e+10, 1.989e+33, 3.826e+33 
        self.sigma, self.c, self.a =  5.67051e-5, 2.99792458e+10, 7.56591e-15
        self.G, self.k_B, self.m_H = 6.67259e-8, 1.380658e-16, 1.673534e-24
        self.pi, self.gamma = 3.141592654e0, 1.6666667e0
        self.tog_bf, self.g_ff = 0.01e0, 1.0e0

        # Rsun = radio del Sol
        # Msun = masa del Sol
        # Lsun = luminosidad del Sol
        # sigma = constante de Stefan-Boltzmann
        # c = velocidad de la luz en el vacío
        # a = 4*sigma/c (constante de presión de radiación)
        # G = constante de gravitación universal
        # K_B = constante de Boltzmann
        # m_H = masa del átomo de hidrógeno
        # pi = 3.141592654
        # gamma = 5/3 (gamma adiabatica para un gas monoatómico)
        # gamrat = gamma/(gamma-1)
        # kPad = P/T**(gamma/(gamma-1)) (constante de un gas adiabático)
        # tog_bf = bound-free opacity constant (ratio of guillotine to gaunt factors)
        # g_ff = free-free opacity gaunt factor (assumed to be unity)



        """ Parametros del modelo """

        # Parametros del modelo
        self.Nstart, self.Nstop, self.igoof, self.ierr = 10, 999, -1, 0
        self.p0, self.t0, self.dlplim = 0.0e0, 0.0e0, 99.9e0 

        # Fracción de masa del ciclo CNO definida a 0.5*Z
        self.XCNO = 0.5*self.Z

        # Conversión de unidades y cálculo del radio solar
        # Ms,Ls,Rs = Masa, luminosidad y radio de la estrella (cgs)
        self.Ms = self.Msolar*self.Msun
        self.Ls = self.Lsolar*self.Lsun
        self.Rs = np.sqrt(self.Ls/(4.e0*self.pi*self.sigma))/(self.Te**2)
        self.Rsolar = self.Rs/self.Rsun

        #  Begin with a very small step size since surface conditions vary
        #  rapidly.
        #
        self.deltar = -self.Rs/1000.0e0
        self.idrflg = 0
        #     
        #  Calculate mean molecular weight mu assuming complete ionization
        #  (see Eq. 10.21).
        #
        self.mu = 1.0e0/(2.0e0*self.X + 0.75e0*self.Y + 0.5e0*self.Z)
        #
        #  Calculate the delimiter between adiabatic convection and radiation
        #  (see Eq. 10.87).
        #
        self.gamrat = self.gamma/(self.gamma - 1.0e0)

    def discretizeVariables(self):

        # Crearemos arreglos de tamaño (999)
        self.size = 999
        N = self.size

        # Creamos los arreglos que vamos a usar para almacenar
        # las soluciones de la estructura estelar
        self.R, self.M_r = np.zeros(N), np.zeros(N)
        self.P, self.T = np.zeros(N), np.zeros(N)
        self.L_r, self.rho = np.zeros(N), np.zeros(N)
        self.kappa, self.epsilon = np.zeros(N), np.zeros(N)

        # Definimos las condiciones de frontera en 
        # la superficie de la estrella
        self.R[N], self.M_r[N] = self.Rs, self.Ms
        self.L_r[N], self.T[N], self.P[N] = self.Ls, self.T0, self.P0

    
    def main(self):

        # FORMATO PARA LA IMPRESIÓN DE MENSAJES

        def formato100():
            formato = "Las fracciones de hidrogeno X y de metales Z no son válidas"
            print(formato)

        def formato200():
            formato = """The variation in mass has become larger than 0.001*Ms leaving the approximation loop before Nstart was reached"""
            print(formato)

        def formato300():
            formato = """The problem occurred in the Runge-Kutta routine"""
            print(formato)

        def formato400(r, rho, M_r, kappa, T, epsilon, P, L_r, Rs, Ms, Ls):
            formato = f"""Values from the previous zone are:
                         r/Rs    = {r/Rs}
                         rho     = {rho} g/cm**3 
                         M_r/Ms  = {M_r/Ms} 
                         kappa   = {kappa} cm**2/g
                         T       = {T} K
                         epsilon = {epsilon} ergs/g/s
                         P       = {P} dynes/cm**2',/,
                         L_r/Ls  = {L_r/Ls}"""
            print(formato)

        def formato1000():
            formato = """A Homogeneous Main-Sequence Model"""
            print(formato)

        def formato2000(Msolar, Mcrat, Rsolar, Rcrat, Lsolar, Lcrat, Te,
            rhocor, X, Tcore, Y, Pcore, Z, epscor, dlpdlt, istop):
            formato = f"""The surface conditions are:          The central conditions are:
            Mtot = {Msolar} Msun            Mc/Mtot     = {Mcrat}
            Rtot = {Rsolar} Rsun            Rc/Rtot     = {Rcrat}
            Ltot = {Lsolar} Lsun            Lc/Ltot     = {Lcrat}
            Teff = {Te} K                   Density     = {rhocor} g/cm**3
            X    = {X}                      Temperature = {Tcore} K
            Y    = {Y}                      Pressure    = {Pcore} dynes/cm**2
            Z    = {Z}                      epsilon     = {epscor} ergs/s/g
                                            dlnP/dlnT   = {dlpdlt}"""
            print(formato)

        def formato2500():
            formato = """"""
            print(formato)

        def formato3000():
            formato = """"""
            print(formato)    

        def formato4000():
            formato = """"""
            print(formato) 

        def formato5000():
            formato = """"""
            print(formato)

        def formato5100():
            formato = """"""
            print(formato)

        def formato5200():
            formato = """"""
            print(formato) 

        def formato5300():
            formato = """"""
            print(formato)

        def formato5400():
            formato = """"""
            print(formato)

        def formato5500():
            formato = """"""
            print(formato)

        def formato5600():
            formato = """"""
            print(formato)

        def formato5700():
            formato = """"""
            print(formato)

        def formato5800():
            formato = """"""
            print(formato)   

        def formato6000():
            formato = """"""
            print(formato)

        def formato7000():
            formato = """"""
            print(formato)

        def formato9000():
            formato = """"""
            print(formato)                                     


        # PROGRAMA PRINCIPAL

        # Validamos las fracciones X y Z
        if(self.Y < 0):
            formato100()

        N = self.size

        if (self.P0 <= 0.0e0 or self.T0 <= 0.0e0):
            self.rho[N] = 0.0e0
            self.kappa[N]  = 0.0e0
            self.epsilon[N] = 0.0e0
        else:
            call EOS(X, Z, XCNO, mu, P(1), T(1), rho(1), kappa(1),
        1        epslon(1), tog_bf, 1 ,ierr)
            if (ierr.ne.0) stop
        end if



    def EOS(self, T, P, rho, kappa, epsilon, tog_bg, izone):

        def formato100(T, P, izone):
            formato = f"""'Something is a little wrong here.
            You are asking me to deal with either a negative temperature'
            or a negative pressure.  I am sorry but that is not in my
            contract! You will have to try again with different initial 
            conditions.
            In case it helps, I detected the problem in zone {izone}
            with the following conditions:
            T = {T} K,
            P = {P} dynes/cm**2 """
            print(formato)

        def formato200(T, P, izone, Prad, Pgas, rho):
              formato = f"""'I am sorry, but a negative density was detected,
                my equation-of-state routine is a bit baffled by this new
                physical system you have created.  The radiation pressure
                is probably too great, implying that the star is unstable.
                Please try something a little less radical next time.
                In case it helps, I detected the problem in zone {izone}
                with the following conditions:
                T       = {T} K ,
                P_total = {P} dynes/cm**2 ,
                P_rad   = {Prad} dynes/cm**2 ,
                P_gas   = {Pgas} dynes/cm**2 ,
                rho     = {rho} g/cm**3 """
              print(formato)


        oneo3, twoo3 = 0.333333333e0, 0.666666667e0

        if (T <= 0.0e0 or P <= 0.0e0):
            self.ierr = 1
            formato100(T, P, izone)

        Prad = self.a*T**4/3.0e0
        Pgas = P - Prad
        rho = (self.mu*self.m_H/self.k_B)*(Pgas/T)

        if (rho < 0.0e0):
            self.ierr = 1
            formato200(T, P, izone, Prad, Pgas, rho)

        # Calculate opacity, including the guillotine-to-gaunt factor ratio;
        # see Novotny (1973), p. 469. k_bf, k_ff, and k_e are the bound-free,
        # free-free, and electron scattering opacities, given by Eqs. (9.19),
        #  (9.20), and (9.21), respectively.

        tog_bf = 2.82e0*(rho*(1.0e0 + self.X))**0.2e0
        k_bf = 4.34e25/tog_bf*self.Z*(1.0e0 + self.X)*rho/T**3.5e0
        k_ff = 3.68e22*self.g_ff*(1.0e0 - self.Z)*(1.0e0 + self.X)*rho/T**3.5e0
        k_e = 0.2e0*(1.0e0 + self.X)
        kappa = k_bf + k_ff + k_e

        # Compute energy generation by the pp chain and the CNO cycle.  These
        # are calculated using Eqs. (10.49) and (10.53), which come from 
        # Fowler, Caughlan, and Zimmerman (1975). The screening factor for the
        # pp chain is calculated as fpp; see Clayton (1968), p. 359ff.

        T6 = T*1.0e-06
        fx = 0.133e0*self.X*np.sqrt((3.0e0 + self.X)*rho)/T6**1.5e0
        fpp = 1.0e0 + fx*self.X
        psipp = 1.0e0 + 1.412e8*(1.0e0/self.X - 1.0e0)*np.exp(-49.98*T6**(-oneo3))
        Cpp = 1.0e0 + 0.0123e0*T6**oneo3 + 0.0109e0*T6**twoo3 + 0.000938e0*T6
        epspp = 2.38e6*rho*self.X*self.X*fpp*psipp*Cpp*T6**(-twoo3)*np.exp(-33.80e0*T6**(-oneo3))
        CCNO = 1.0e0 + 0.0027e0*T6**oneo3 - 0.00778e0*T6**twoo3 - 0.000149e0*T6
        epsCNO = 8.67e27*rho*self.X*self.XCNO*CCNO*T6**(-twoo3)*np.exp(-152.28e0*T6**(-oneo3))
        epsilon = epspp + epsCNO
 
