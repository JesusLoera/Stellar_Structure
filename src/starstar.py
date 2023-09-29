
""" Librerías """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

""" Clase Stellar Structure"""

class StellarStructure():

    def __init__(self, Msolar, Lsolar, Te, X, Z, nsh=999):

 
        """ Constrcutor con los parámetros de la estrella """

        # Características de la estrella
        self.Msolar = Msolar      # Masa (unidades solares)
        self.Lsolar = Lsolar      # Luminosidad (unidades solares)
        self.Te = Te              # Temperatura efectiva (kelvins)
        self.X = X                # Fracción de hidrogeno (adimensional)
        self.Z = Z                # Fracción de metales (adimensional)
        self.Y = 1.0e0 - X - Z    # Fracción de helio (adimensional)
        self.nsh = nsh

        # Validamos las fracciones X y Z
        while((self.Y < 0)):
            print("Las fracciones de hidrogeno X y de metales Z no son válidas")
            self.X = float(input(' Enter the mass fraction of hydrogen (X): '))
            self.Z = float(input(' Enter the mass fraction of metals (Z): '))
            self.Y = 1.e0 - self.X - self.Z

       

    def initializeVariables(self):

        # Fracción de masa del ciclo CNO definida a 0.5*Z
        self.XCNO = 0.5*self.Z

        #  Calculate mean molecular weight mu assuming complete ionization
        #  (see Eq. 10.21).
        self.mu = 1.0e0/(2.0*self.X + 0.75*self.Y + 0.5*self.Z)

        #  Calculate the delimiter between adiabatic convection and radiation
        #  (see Eq. 10.87).
        self.gamrat = self.gamma/(self.gamma - 1.0e0)

        # Crearemos arreglos de tamaño (999)
        N = self.nsh

        # Creamos los arreglos que vamos a usar para almacenar
        # las soluciones de la estructura estelar
        self.r, self.P = np.zeros(N, float), np.zeros(N, float)
        self.M_r, self.L_r = np.zeros(N, float), np.zeros(N, float)
        self.T, self.rho = np.zeros(N, float), np.zeros(N, float)
        self.kappa, self.epsilon = np.zeros(N, float), np.zeros(N, float)
        self.tog_bf, self.dlPdlT = np.zeros(N, float), np.zeros(N, float)

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

        # initialize variables used for 4 structure equation derivatives

        self.f_im1=np.zeros(4,float)
        self.dfdr=np.zeros(4,float)
        self.f_i=np.zeros(4,float)


    """ Método EOS """
    # Cálcula los valores de la densidad, la opacidad, el 
    # guillotine-to-gaunt factor ratio y la razón de generación 
    # de la energía en un radio r. 

    def EOS(self, P, T, izone):

        def formato100(P, T, izone):
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

        def formato200(P, T, izone, Prad, Pgas, rho):
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

        if ((T <= 0.0e0) or (P <= 0.0e0)):
            formato100(P, T, izone)
            return (0.0, 0.0, 0.0, 0.0, 1)


        Prad = self.a*T**4/3.0e0
        Pgas = P - Prad
        rho = (self.mu*self.m_H/self.k_B)*(Pgas/T)

        if (rho < 0.0e0):
            formato200(P, T, izone, Prad, Pgas, rho)
            return (rho, 0.0 , 0.0 ,0.0, 1)

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

        return (rho, kappa, epsilon, tog_bf, 0)
    
    def STARTMDL(self, r_i, M_ri, L_ri, tog_bf, irc):

        r = r_i + self.deltar
        M_rip1 = M_ri
        L_rip1 = L_ri

        #  This is the radiative approximation (neglect radiation pressure
        #  and electron scattering opacity)# see Prialnik Eq. 5.1, 5.3 and Sec. 3.7 or C&O Eqs. (H.1), (H.2), (9.19),
        #  and (9.20).

        if (irc == 0):
            T_ip1 = self.G*M_rip1*self.mu*self.m_H/(4.25e0*self.k_B)*(1.0e0/r - 1.0e0/self.Rs)
            A_bf = 4.34e25*self.Z*(1.0e0 + self.X)/tog_bf
            A_ff = 3.68e22*self.g_ff*(1.0e0 - self.Z)*(1.0e0 + self.X)
            Afac = A_bf + A_ff
            P_ip1 = np.sqrt((1.0e0/4.25e0)*(16.0e0/3.0e0*np.pi*self.a*self.c)*(self.G*M_rip1/L_rip1)*(self.k_B/(Afac*self.mu*self.m_H)))*T_ip1**4.25e0
        #  This is the convective approximation# see Prialnik Sec 6.5, 6.6 or C&O Eqs. (H.3) and (10.75).
        else:
            T_ip1 = self.G*M_rip1*self.mu*self.m_H/self.k_B*(1.0e0/r - 1.0e0/self.Rs)/self.gamrat
            P_ip1 = self.kPad*T_ip1**self.gamrat

        return r,P_ip1, M_rip1, L_rip1, T_ip1



 

    
    def main(self):

        # PROGRAMA PRINCIPAL

         # FORMATO PARA LA IMPRESIÓN DE MENSAJES

        def formato200():
            formato = """The variation in mass has become larger than 0.001*Ms leaving the approximation loop before Nstart was reached"""
            print(formato)

        def formato300():
            formato = """The problem occurred in the Runge-Kutta routine"""
            print(formato)

        def formato400(r, rho, M_r, kappa, T, epsilon, P, L_r):
            formato = f"""Values from the previous zone are:
                         r/Rs    = {r/self.Rs}
                         rho     = {rho} g/cm**3 
                         M_r/Ms  = {M_r/self.Ms} 
                         kappa   = {kappa} cm**2/g
                         T       = {T} K
                         epsilon = {epsilon} ergs/g/s
                         P       = {P} dynes/cm**2',/,
                         L_r/Ls  = {L_r/self.Ls}"""
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

        def formato2500(Ms):
            formato = f"""Notes:
            (1) Mass is listed as Qm = 1.0 - M_r/Mtot, where Mtot = {Ms} g
            (2) Convective zones are indicated by c, radiative zones by r
            (3) dlnP/dlnT may be limited to +99.9 or -99.9; if so it is labeled by *)"""
            print(formato)

        def formato4000():
            formato = """Sorry to be the bearer of bad news, but... Your model has some problems"""
            print(formato) 

        def formato5000():
            formato = """The number of allowed shells has been exceeded"""
            print(formato)

        def formato5100():
            formato = """"""
            print(formato)

        def formato5200(rho, istop):
            formato = f"""The core density seems a bit off,density should increase smoothly toward the center.The density of the last zone calculated was rho = {rho} gm/cm**3"""
            print(formato) 

        def formato5300():
            formato = """It looks like you will need a degenerate neutron gas and general relativity to solve this core.  Who do you think I am, Einstein?"""
            print(formato)

        def formato5400(epsilon, istop):
            formato = f"""The core epsilon seems a bit off epsilon should vary smoothly near the center. The value calculated for the last zone was eps = {epsilon} ergs/g/s"""
            print(formato)

        def formato5500(T, istop):
            formato = f"""Your extrapolated central temperature is too low a little more fine tuning ought to do it.
            The value calculated for the last zone was T = {T} K"""
            print(formato)

        def formato5600():
            formato = """You created a star with a hole in the center!"""
            print(formato)

        def formato5700():
            formato = """This star has a negative central luminosity!"""
            print(formato) 

        def formato5800():
            formato = """You hit the center before the mass and/or luminosity were depleted!"""
            print(formato)

        def formato6000():
            formato = """It looks like you are getting close, however, there are still a few minor errors"""
            print(formato)

        def formato7000():
            formato = """CONGRATULATIONS, I THINK YOU FOUND IT! However, be sure to look at your model carefully."""
            print(formato)

        def formato9000():
            formato = """***** The integration has been completed *****
            The model has been stored in starmodl.dat'"""
            print(formato)  

        # Definicón del constantes
        self.Rsun, self.Msun, self.Lsun = 6.9599e+10, 1.989e+33, 3.826e+33 
        self.sigma, self.c, self.a =  5.67051e-5, 2.99792458e+10, 7.56591e-15
        self.G, self.k_B, self.m_H = 6.67259e-8, 1.380658e-16, 1.673534e-24
        self.pi, self.gamma = 3.141592654e0, 1.6666667e0
        self.tog_bf0, self.g_ff = 0.01, 1.0e0

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


        # Parametros del modelo
        self.Nstart, self.Nstop, self.igoof, self.ierr = 10, self.nsh, -1, 0
        self.P0, self.T0, self.dlPlim, self.debug = 0.0, 0.0, 99.9, 0 


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

        self.initializeVariables()
        initsh = 0

        # Definimos las condiciones de frontera en 
        # la superficie de la estrella
        self.r[initsh], self.M_r[initsh] = self.Rs, self.Ms
        self.L_r[initsh], self.T[initsh] = self.Ls, self.T0
        self.P[initsh], self.tog_bf[initsh] = self.P0, self.tog_bf0

        if (self.P0 <= 0.0) or (self.T0 <= 0.0):
            self.rho[initsh]    = 0.0
            self.kappa[initsh]  = 0.0
            self.epslon[initsh] = 0.0
            self.tog_bf[initsh] = 0.01
        else:
            self.rho[initsh], self.kappa[initsh], self.epsilon[initsh], self.tog_bf[N], self.ierr = self.EOS(self.P[initsh], self.T[initsh], initsh)
            if (self.ierr != 0):
                print ("we're stopping now")
                istop=0

        #  Apply approximate surface solutions to begin the integration,
        #  assuming radiation transport in the outermost zone (the do 20 loop).
        #  irc = 0 for radiation, irc = 1 for convection.
        #  Assume arbitrary initial values for kPad, and dlPdlT.
        #  dlPdlT = dlnP/dlnT (see Prialnik Eq. 6.28 or C&O Eq. 10.87)

        self.kPad = 0.3e0
        self.irc = 0
        self.dlPdlT[initsh] = 4.25e0
        for i in range(0, self.Nstart):

            ip1 = i + 1

            self.r[ip1], self.P[ip1], self.M_r[ip1], self.L_r[ip1], self.T[ip1]= self.STARTMDL(self.r[i], self.M_r[i], self.L_r[i], self.tog_bf[i], ip1)

            self.rho[ip1], self.kappa[ip1], self.epsilon[ip1], self.tog_bf[ip1],self.ierr = self.EOS(self.P[ip1], self.T[ip1], ip1)


            if (self.ierr != 0):
                formato400(self.r[i], self.rho[i], self.M_r[i], 
                           self.kappa[i], self.T[i], self.epsilon[i],
                           self.P[i],
                           self.L_r[i])      
                break                       

            
            #  Determine whether convection will be operating in the next zone
            #  calculating dlnP/dlnT numerically between zones i and i+1 [ip1].
            #  Update the adiabatic gas constant if necessary. 

            if (i > initsh):
                self.dlPdlT[ip1] = np.log(self.P[ip1]/self.P[i])/np.log(self.T[ip1]/self.T[i])
            else:
                self.dlPdlT[ip1] = self.dlPdlT[i]

            if (self.dlPdlT[ip1] < self.gamrat):
                irc = 1
            else:
                irc = 0
                self.kPad = self.P[ip1]/self.T[ip1]**self.gamrat

        # PROGRAMA PRINCIPAL


        # if (self.P0 <= 0.0e0 or self.T0 <= 0.0e0):
        #     self.rho[N] = 0.0e0
        #     self.kappa[N]  = 0.0e0
        #     self.epsilon[N] = 0.0e0
        # else:
            
        #     self.rho[N], self.kappa[N], self.epsilon[N], self.tog_bf[N], ierr = self.EOS(X, Z, XCNO, mu, P[initsh], T[initsh], 0 ,cst)self.EOS(self.P[N], self.T[N], self.rho[N], self.kappa[N],
        #              self.epsilon[N], self.tog_bf, 1 , self.ierr)
        #     if (self.ierr != 0):



star = StellarStructure(1, 1, 6000, 0.9,  18)