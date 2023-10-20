"""
Se ha reescrito el código StatStar_python3.py a un nuevo código 
empleando clases, con la finalidad de que el código sea un poco
más oordenado y más limpio de usar cuando se quiere realizar 
una gran cantidad de cálculos de diferentes estructuras estelares.

Por: Jesús Loera

"""

# Encontré un poco confuso el código orginal pero al final pudo 
# correr y replicar el código original, pienso seguir trabajando
# en la limpieza de este código porque considero que tiene aún
# un estilo muy "Fortran", cuando podría escribirse de una
# manera más "Pythonica".


""" Librerías """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

""" Clase Stellar Structure"""

class StellarStructure():

    def __init__(self, Msolar, Lsolar, Te, X, Z, filename, nsh=999):

 
        """ Constructor con los parámetros de la estrella """

        # Características de la estrella
        self.Msolar = Msolar      # Masa (unidades solares)
        self.Lsolar = Lsolar      # Luminosidad (unidades solares)
        self.Te = Te              # Temperatura efectiva (kelvins)
        self.X = X                # Fracción de hidrogeno (adimensional)
        self.Z = Z                # Fracción de metales (adimensional)
        self.Y = 1.0e0 - X - Z    # Fracción de helio (adimensional)
        self.nsh = nsh
        self.filename = filename

        # Validamos las fracciones X y Z
        while((self.Y < 0)):
            print("Las fracciones de hidrogeno X y de metales Z no son válidas")
            self.X = float(input(' Enter the mass fraction of hydrogen (X): '))
            self.Z = float(input(' Enter the mass fraction of metals (Z): '))
            self.Y = 1.e0 - self.X - self.Z

       

    def initializeVariables(self):

        # Fracción de masa del ciclo CNO definida a 0.5*Z
        self.XCNO = self.Z/2.0e0

        #  Calcula el peso molecular medio mu asumiendo una ionización completa
        #  (vea Eq. 10.21).
        self.mu = 1.0e0/(2.0*self.X + 0.75*self.Y + 0.5*self.Z)

        #  Calcula el delimitador entre convección adiabatica y radiación
        #  (vea Eq. 10.87).
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

        #inicializar variables utilizadas para 4 derivadas de ecuaciones de estructura

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


        oneo3=0.333333333e0
        twoo3=0.666666667e0

        if ((T < 0.0e0) or (P < 0.0e0)):
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

        #  Esta es la aproximación radiativa (neglect radiation pressure
        #  and electron scattering opacity)# see Prialnik Eq. 5.1, 5.3 and Sec. 3.7 or C&O Eqs. (H.1), (H.2), (9.19),
        #  and (9.20).

        if (irc == 0):
            T_ip1 = self.G*M_rip1*self.mu*self.m_H/(4.25e0*self.k_B)*(1.0e0/r - 1.0e0/self.Rs)
            A_bf = 4.34e25*self.Z*(1.0e0 + self.X)/tog_bf
            A_ff = 3.68e22*self.g_ff*(1.0e0 - self.Z)*(1.0e0 + self.X)
            Afac = A_bf + A_ff
            P_ip1 = np.sqrt((1.0e0/4.25e0)*(16.0e0/3.0e0*np.pi*self.a*self.c)*(self.G*M_rip1/L_rip1)*(self.k_B/(Afac*self.mu*self.m_H)))*T_ip1**4.25e0
        # Esta es la aproximación convectiva
        # see Prialnik Sec 6.5, 6.6 or C&O Eqs. (H.3) and (10.75).
        else:
            T_ip1 = self.G*M_rip1*self.mu*self.m_H/self.k_B*(1.0e0/r - 1.0e0/self.Rs)/self.gamrat
            P_ip1 = self.kPad*T_ip1**self.gamrat

        return r,P_ip1, M_rip1, L_rip1, T_ip1

    #  Los siguientes cuatro métodos calculan los gradientes de presión, masa, luminosidad y temperatura en un radio r.

    def dPdr(self, r, M_r, rho):
        return -self.G*rho*M_r/r**2
    
    def dMdr(self, r, rho):
        return (4.0e0*np.pi*rho*r**2)
    
    def dLdr(self, r, rho, epsilon):
      return (4.0e0*np.pi*rho*epsilon*r**2)

    def dTdr(self, r, M_r, L_r, T, rho, kappa, irc):
        if (irc == 0):
            return (-(3.0e0/(16.0e0*np.pi*self.a*self.c))*kappa*rho/T**3*L_r/r**2)
        # Este es el gradiente de temperatura convectiva adiabatica (Prialnik Eq. 6.29 or C&O Eq. 10.81).
        else:
            return (-1.0e0/self.gamrat*self.G*M_r/r**2*self.mu*self.m_H/self.k_B)
        
    # Método FUNDEQ(r, f, irc, izone)

    def FUNDEQ(self, r, f, irc, izone):

        dfdr=np.zeros(4)
        P   = f[0]
        M_r = f[1]
        L_r = f[2]
        T   = f[3]
        rho, kappa, epsilon, tog_bf, ierr = self.EOS(P, T, izone)
        dfdr[0] = self.dPdr(r, M_r, rho)
        dfdr[1] = self.dMdr(r, rho)
        dfdr[2] = self.dLdr(r, rho, epsilon)
        dfdr[3] = self.dTdr(r, M_r, L_r, T, rho, kappa, irc)
        return (dfdr,ierr)


    #
    # Algoritmo de Runge-kutta
    #
    def RUNGE(self, f_im1, dfdr, r_im1, deltar, irc, izone):

        f_temp=np.zeros(4)
        f_i=np.zeros(4)

        dr12 = deltar/2.0e0
        dr16 = deltar/6.0e0
        r12  = r_im1 + dr12
        r_i  = r_im1 + deltar

        # Calcular derivadas intermedias a partir de las cuatro estelares fundamentales encontradas en el método FUNDEQ.

        for i in range(0,4):
            f_temp[i] = f_im1[i] + dr12*dfdr[i]

        df1, ierr = self.FUNDEQ(r12, f_temp, irc, izone)
        if (ierr != 0):
            return f_i,ierr

        for i in range(0,4):
            f_temp[i] = f_im1[i] + dr12*df1[i]

        df2, ierr = self.FUNDEQ(r12, f_temp, irc, izone)

        if (ierr != 0):
            return f_i,ierr

        for i in range(0,4):
            f_temp[i] = f_im1[i] + deltar*df2[i]

        df3, ierr=self.FUNDEQ(r_i, f_temp, irc, izone)
        if (ierr != 0):
            return f_i,ierr

        # Calcula las variables en la siguiente capa (i + 1).

        for i in range(0,4):
            f_i[i] = f_im1[i] + dr16*(dfdr[i] + 2.0e0*df1[i] + 2.0e0*df2[i] + df3[i])

        return f_i,0


    # Programa principal
    
    def main(self):
         
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

        # Definición de constantes
        self.Rsun, self.Msun, self.Lsun = 6.9599e10, 1.989e33, 3.826e33 
        self.sigma, self.c, self.a =  5.67051e-5, 2.99792458e10, 7.56591e-15
        self.G, self.k_B, self.m_H = 6.67259e-8, 1.380658e-16, 1.673534e-24
        self.pi, self.gamma = 3.141592654e0, 5.0e0/3
        self.tog_bf0, self.g_ff = 0.01, 1.0e0

        # deltar = paso de integración para el radio
        # idrflg = flag que define el tamaño del paso de integración
        #        = 0 (tamaño de paso inical para la superficie de Rs/1000.)
        #        = 1 (tamaño de paso estándar de Rs/100.)
        #        = 2 (tamaño de paso para el núcleo de Rs/5000.)
        #  
        # Nstart = número de pasos para los que deben utilizarse ecuaciones de partida
        #          (Se supone que la zona más externa es radiativa)
        # Nstop = Número máximo de zonas permitidas en la estrella.
        # Igoof = flag con la condición final del modelo
        #       = -1 (número de zonas excedido; también el valor inicial)
        #       =  0 (buen modelo)
        #       =  1 (la densidad del núcleo es extrema)
        #       =  2 (la luminosidad del núcleo es extrema)
        #       =  3 (la temperatura extrapolada del núcleo es muy baja)
        #       =  4 (la masa se volvió negativa antes de alcanzar el núcleo)
        #       =  5 (la luminosidad se volció negativa antes de alcanzar el núcleo)
        # T0, P0 = temperatura y presión en la superficie (se asume T0 = P0 = 0)


        # Parametros del modelo
        self.Nstart, self.Nstop, self.Igoof, self.ierr = 10, self.nsh, -1, 0
        self.P0, self.T0, self.dlPlim, self.debug = 0.0, 0.0, 99.9, 0 


        # Conversión de unidades y cálculo del radio solar
        # Ms,Ls,Rs = Masa, luminosidad y radio de la estrella (cgs)
        self.Ms = self.Msolar*self.Msun
        self.Ls = self.Lsolar*self.Lsun
        self.Rs = np.sqrt(self.Ls/(4.e0*self.pi*self.sigma))/(self.Te**2)
        self.Rsolar = self.Rs/self.Rsun

        #  Empezamos con un muy pequeño paso de integración débido a que las condiciones en la superficie varían muy rápido

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
            self.epsilon[initsh] = 0.0
            self.tog_bf[initsh] = 0.01
        else:
            self.rho[initsh], self.kappa[initsh], self.epsilon[initsh], self.tog_bf[initsh], self.ierr = self.EOS(self.P[initsh], self.T[initsh], initsh)
            if (self.ierr != 0):
                print ("we're stopping now")
                istop=0

        # Aplicar soluciones superficiales aproximadas para comenzar la integración,
        # suponiendo transporte de radiación en la zona más externa (el bucle do 20).
        # irc = 0 para radiación, irc = 1 para convección.
        # Supongamos valores iniciales arbitrarios para kPad y dlPdlT.
        # dlPdlT = dlnP/dlnT (consulte la ecuación de Prialnik 6.28 o la ecuación de C&O 10.87)

        self.kPad = 0.3e0
        self.irc = 0
        self.dlPdlT[initsh] = 4.25e0
        for i in range(0, self.Nstart):

            ip1 = i + 1

            self.r[ip1], self.P[ip1], self.M_r[ip1], self.L_r[ip1], self.T[ip1]= self.STARTMDL(self.r[i], self.M_r[i], self.L_r[i], self.tog_bf[i], self.irc)

            self.rho[ip1], self.kappa[ip1], self.epsilon[ip1], self.tog_bf[ip1],self.ierr = self.EOS(self.P[ip1], self.T[ip1], ip1)


            if (self.ierr != 0):
                formato400(self.r[i], self.rho[i], self.M_r[i], 
                           self.kappa[i], self.T[i], self.epsilon[i],
                           self.P[i],
                           self.L_r[i])      
                break                       

            
            # Determinar si la convección funcionará en la siguiente zona
            # calculando dlnP/dlnT numéricamente entre las zonas i e i+1 [ip1].
            # Actualice la constante adiabática del gas si es necesario.

            if (i > initsh):
                self.dlPdlT[ip1] = np.log(self.P[ip1]/self.P[i])/np.log(self.T[ip1]/self.T[i])
            else:
                self.dlPdlT[ip1] = self.dlPdlT[i]

            if (self.dlPdlT[ip1] < self.gamrat):
                irc = 1
            else:
                irc = 0
                self.kPad = self.P[ip1]/self.T[ip1]**self.gamrat

            # Pruebe para ver si la suposición masa constante en la superficie sigue siendo válido.

            deltaM = self.deltar*self.dMdr(self.r[ip1], self.rho[ip1])
            self.M_r[ip1] = self.M_r[i] + deltaM
            if (np.abs(deltaM) > (0.001e0*self.Ms)):
                if (ip1 > 1):
                    ip1 = ip1 - 1
                    print(' The variation in mass has become larger than 0.001*Ms')
                    print(' leaving the approximation loop before Nstart was reached')
                    break

            
        Nsrtp1 = ip1 + 1

        if (self.ierr != 0):    
            # salir si hemos llegado a este punto después de un error en la inicialización
            self.Nstop=Nsrtp1-1
            istop=self.Nstop

        for i in range(Nsrtp1, self.Nstop):
            im1 = i - 1
            self.f_im1[0] = self.P[im1]
            self.f_im1[1] = self.M_r[im1]
            self.f_im1[2] = self.L_r[im1]
            self.f_im1[3] = self.T[im1]
            self.dfdr[0]  = self.dPdr(self.r[im1], self.M_r[im1], self.rho[im1])
            self.dfdr[1] = self.dMdr(self.r[im1], self.rho[im1])
            self.dfdr[2] = self.dLdr(self.r[im1], self.rho[im1], self.epsilon[im1])
            self.dfdr[3] = self.dTdr(self.r[im1], self.M_r[im1], self.L_r[im1], self.T[im1], self.rho[im1], self.kappa[im1], irc)
            f_i,ierr=self.RUNGE(self.f_im1, self.dfdr, self.r[im1], self.deltar, irc, i)

            if (ierr != 0):
                formato300()
                formato400(self.r[im1], self.rho[im1], self.M_r[im1], self.kappa[im1], self.T[im1], self.epsilon[im1], self.P[im1], self.L_r[im1])
                break

            # Actualizar los parámetros estelares para la siguiente zona, incluyendo la adición de dr al radio antiguo (nótese que dr < 0 ya que la integración es hacia el interior).

            self.r[i]   = self.r[im1] + self.deltar
            self.P[i]   = f_i[0]
            self.M_r[i] = f_i[1]
            self.L_r[i] = f_i[2]
            self.T[i]   = f_i[3]

            # Calcule la densidad, opacidad y tasa de generación de energía para esta zona.

            self.rho[i], self.kappa[i], self.epsilon[i], self.tog_bf[i], self.ierr = self.EOS(self.P[i], self.T[i], i)

            if (ierr != 0):
                formato400(self.r[im1], self.rho[im1], self.M_r[im1], self.kappa[im1], self.T[im1], self.epsilon[im1], self.P[im1], self.L_r[im1])
                istop = i
                break

            if (self.debug == 1): 
                print (i, self.r[i], self.M_r[i], self.L_r[i], self.T[i], self.P[i], self.rho[i], self.kappa[i], self.epsilon[i], self.tog_bf[i])

            # Determinar si la convección funcionará en la siguiente zona
            # calculando dlnP/dlnT y comparándolo con gamma/(gamma-1)
            # (consulte la ecuación de Prialnik 6.28 o la ecuación de C&O 10.87). Configure la bandera de convección apropiadamente.

            self.dlPdlT[i] = np.log(self.P[i]/self.P[im1])/np.log(self.T[i]/self.T[im1])

            if (self.dlPdlT[i] < self.gamrat):
                irc = 1
            else:
                irc = 0

            # Comprueba si se ha alcanzado el centro.  Si es así, establece Igoof y
            # estimar las condiciones centrales rhocor, epscor, Pcore, y Tcore.
            # La densidad central se estima que es la densidad media de la
            # bola central restante, la presión central se determina
            # utilizando la expansión de Taylor en el centro (Prialnik - Ejercicio. 5.1; CO Ec. H.4)
            # y el valor central de la tasa de generación
            # de generación de energía se calcula como la luminosidad interior restante
            # luminosidad interior dividida por la masa de la bola central.  Finalmente, la # temperatura central se calcula
            # temperatura central se calcula aplicando la ley de los gases ideales
            # (despreciando la presión de radiación).

            if ((self.r[i] <= np.abs(self.deltar)) and ((self.L_r[i] >= (0.1e0*self.Ls)) or (self.M_r[i] >= (0.01e0*self.Ms)))):
                # Llega al centro antes de que se agote la masa/luminosidad
                self.Igoof = 6
            elif (self.L_r[i] <= 0.0e0):
                # Obtenido luminosidad central negativa
                self.Igoof = 5
                # rho: expansión de taylor en el centro
                rhocor = self.M_r[i]/(4.0e0/3.0e0*np.pi*self.r[i]**3)
                if (self.M_r[i] != 0.0e0):
                    # razón de generación de energía
                    epscor = self.L_r[i]/self.M_r[i]
                else:
                    epscor = 0.0e0
                    # P:  Expanxión de Taylor en el centro
                    Pcore = self.P[i] + 2.0e0/3.0e0*np.pi*self.G*rhocor**2*self.r[i]**2
                    # Se asume un gas ideal
                    Tcore = Pcore*self.mu*self.m_H/(rhocor*self.k_B)
            elif (self.M_r[i] <= 0.0e0):
                self.Igoof  = 4  # El módelo tiene un hoyo en el centro (densidad negativa!)
                rhocor = 0.0e0
                epscor = 0.0e0
                Pcore  = 0.0e0
                Tcore  = 0.0e0
            elif ((self.r[i] < (0.02e0*self.Rs)) and ((self.M_r[i] < (0.01e0*self.Ms)) and ((self.L_r[i] < 0.1e0*self.Ls)))):
                #  if we've reached <2% star's radius,
                #<1% mass enclosed and <10% luminosity then....
                # rho: Taylor expansion at center
                rhocor = self.M_r[i]/(4./3.*np.pi*self.r[i]**3)
                # set maximum reasonable core mass
                rhomax = 10.0e0*(self.rho[i]/self.rho[im1])*self.rho[i]
                epscor = self.L_r[i]/self.M_r[i]
                # P: Taylor expansion at center
                Pcore  = self.P[i] + 2.0e0/3.0e0*np.pi*self.G*rhocor**2*self.r[i]**2
                # Assume ideal gas
                Tcore  = Pcore*self.mu*self.m_H/(rhocor*self.k_B)
                # In general, these should all produce values
                # that rise towards center (but not too high)
                if ((rhocor < self.rho[i]) or (rhocor > rhomax)):
                    # rho is off either large or small
                    self.Igoof = 1
                elif (epscor < self.epsilon[i]):
                    # energy generation rate a bit off (low)
                    self.Igoof = 2
                elif (Tcore < self.T[i]):
                    # Temperature a bit off (low)
                    self.Igoof = 3
                else:
                    # number of allowed shells has been exceeded
                    self.Igoof = 0

            if (self.Igoof != -1):
                istop = i
                break

            #  Is it time to change the step size?

            if ((self.idrflg == 0) and (self.M_r[i] < (0.99e0*self.Ms))):
               self.deltar = (-1.0)*self.Rs/100.0e0
               self.idrflg = 1

            if ((self.idrflg == 1) and (self.deltar >= (0.5*self.r[i]))):
                self.deltar = (-1.0)*self.Rs/5000.0e0
                self.idrflg = 2

            istop = i

        #  Generate warning messages for the central conditions.

        rhocor = self.M_r[istop]/(4.0e0/3.0e0*np.pi*self.r[istop]**3)
        epscor = self.L_r[istop]/self.M_r[istop]
        Pcore  = self.P[istop] + 2.0e0/3.0e0*np.pi*self.G*rhocor**2*self.r[istop]**2
        Tcore  = Pcore*self.mu*self.m_H/(rhocor*self.k_B)

        if  (self.Igoof != 0):
            if (self.Igoof == -1):
                formato4000()
                formato5000()

            if (self.Igoof == 1):
                formato6000()
                formato5200(self.rho[istop], istop)
                print (rhocor,rhomax)

            if (rhocor > 1e10):
                formato5300()

            if (self.Igoof == 2):
                formato6000()
                formato5400(self.epsilon[istop], istop)

            if (self.Igoof == 3):
                formato6000()
                formato5500(self.T[istop], istop)

            if (self.Igoof == 4):
                formato4000()
                formato5600()

            if (self.Igoof == 5):
                formato4000()
                formato5700()

            if (self.Igoof == 6):
                formato4000()
                formato5800()
        else:
            formato7000()

        #  Print the central conditions.  If necessary, set limits for the
        #  central radius, mass, and luminosity if necessary, to avoid format
        #  field overflows.

        Rcrat = self.r[istop]/self.Rs
        if (Rcrat < -9.999e0): 
            Rcrat = -9.999e0
        
        Mcrat = self.M_r[istop]/self.Ms
        if (Mcrat < -9.999e0): 
            Mcrat = -9.999e0
        
        Lcrat = self.L_r[istop]/self.Ls
        if (Lcrat < -9.999e0): 
            Lcrat = -9.999e0

        f=open(self.filename, 'w')

        f.write('A Homogeneous Main-Sequence Model\n')
        f.write(' The surface conditions are:        The central conditions are:\n')
        f.write(' Mtot = {0:13.6E} Msun          Mc/Mtot     = {1:12.5E}\n'.format(self.Msolar,Mcrat))
        f.write(' Rtot = {0:13.6E} Rsun          Rc/Rtot     = {1:12.5E}\n'.format(self.Rsolar,Rcrat))
        f.write(' Ltot = {0:13.6E} Lsun          Lc/Ltot     = {1:12.5E}\n'.format(self.Lsolar,Lcrat))
        f.write(' Teff = {0:13.6E} K             Density     = {1:12.5E}\n'.format(self.Te,rhocor))
        f.write(' X    = {0:13.6E}               Temperature = {1:12.5E}\n'.format(self.X,Tcore))
        f.write(' Y    = {0:13.6E}               Pressure    = {1:12.5E} dynes/cm**2\n'.format(self.Y,Pcore))
        f.write(' Z    = {0:13.6E}               epsilon     = {1:12.5E} ergs/s/g\n'.format(self.Z,epscor))
        f.write('                                    dlnP/dlnT   = {0:12.5E}\n'.format(self.dlPdlT[istop]))

        f.write('Notes:\n')
        f.write(' (1) Mass is listed as Qm = 1.0 - M_r/Mtot, where Mtot = {0:13.6}\n'.format(self.Msun))
        f.write(' (2) Convective zones are indicated by c, radiative zones by r\n')
        f.write(' (3) dlnP/dlnT may be limited to +99.9 or -99.9# if so it is\n')
        f.write(' labeled by *\n')

        #  Print data from the center of the star outward, labeling convective
        #   or radiative zones by c or r, respectively.  If abs(dlnP/dlnT)
        #  exceeds 99.9, set a print warning flag (*) and set the output limit
        #  to +99.9 or -99.9 as appropriate to avoid format field overflows.

        f.write('   r        Qm       L_r       T        P        rho      kap      eps     dlPdlT\n')

        for ic in range(0,istop+1):
            i = istop - ic
            Qm = 1.0e0 - self.M_r[i]/self.Ms    # Total mass fraction down to radius


            if (self.dlPdlT[i] < self.gamrat):
                rcf = 'c'
            else:
                rcf = 'r'
            if (np.abs(self.dlPdlT[i]) > self.dlPlim):
                self.dlPdlT[i] = np.copysign(self.dlPlim, self.dlPdlT[i])
                clim = '*'
            else:
                clim = ' '
            s='{0:7.2E} {1:7.2E} {2:7.2E} {3:7.2E} {4:7.2E} {5:7.2E} {6:7.2E} {7:6.2E}{8:1s}{9:1s} {10:5.1f}\n'.format(self.r[i], Qm, self.L_r[i], self.T[i], self.P[i], self.rho[i], self.kappa[i], self.epsilon[i], clim, rcf, self.dlPdlT[i])
            f.write(s)

        # Output to screen

        formato9000()

        return self.Igoof,ierr,istop
    
Msolar = 0.75
Lsolar = 0.189
Te = 3788.5
X = 0.7
Z = 0.008
filename = "star1.dat"
star1 = StellarStructure(Msolar, Lsolar, Te, X, Z, filename)
star1.main()