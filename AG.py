# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 23:00:33 2024

@author: Pablo Ramos
"""
################################################################################################################
### IMPORTACIONES
################################################################################################################
import cantera as ct 
import numpy as np 
from matplotlib import pyplot as plt 
import time as tm
from scipy.integrate import trapz 
from matplotlib.ticker import MultipleLocator
import pandas as pd
import random
import csv
import datetime


especies_inyeccion= ['CH3OH', 'C2H5OH', 'C3H5OH', 'PC4H8OH', 'C2H5O', 'CH3OCH3']
especies_admision=['NC7H16']
historial_origen = [] 




################################################################################################################
### DEFINICIÓN FUNCIONES A EMPLEAR
################################################################################################################


def creacion_individuos(num_individuos, longitud, especies):
    """
    Función para la creación de nuevos individuos. 
    -num_individuos indica el tamaño de la población a crear 
    -longitud indica el tamaño de la cadena binaria para cada uno de los individuos
    
    La variable longitud vendrá determinada principalmente por las especies contempladas 
    """
    individuos_unicos = set()
    nuevos_individuos = []
    
    # Mientras no existan suficientes individuos únicos, se sigue generando nuevos individuos
    while len(individuos_unicos) < num_individuos:
        nuevo_individuo = tuple(np.random.randint(2, size=longitud-2*especies))  # Generar un nuevo individuo como una tupla de bits aleatorios
        presencia_especies = nuevo_individuo[15:] 
        
        # Calcula la suma de los bits de presencia/ausencia
        suma = sum(presencia_especies)
        
        # Verifica que el número de especies presentes sea menor o igual a 3 
        if 1 <= suma <= 3:
            #Asigna los bits de composición según los bits de presencia/ausencia
            for i, presencia in enumerate(presencia_especies):
                if presencia == 1: 
                    composicion = tuple(np.random.randint(2, size=2))
                    nuevo_individuo += composicion
                else:
                    composicion = (0, 0)
                    nuevo_individuo += composicion
                    
            individuos_unicos.add(nuevo_individuo) # Agrega el nuevo individuo al conjunto de individuos únicos
            
    
    # Convertir el conjunto de individuos únicos de nuevo a una lista de listas de bits
    for individuo in individuos_unicos:
        nuevos_individuos.append(list(individuo))
        
    return nuevos_individuos


"""
Funciones para segmentar la cadena de bits del individuo en función de la variable a la que representan 
"""
def obtener_binario_dosado(individuo):
    genotipo_dosado = individuo[0:4]    #Extrae los primero 4 bits dedicados a la codificación del dosado
    return genotipo_dosado
def obtener_binario_RPM(individuo):
    genotipo_RPM = individuo[4:8]       #Extrae los 4 bits dedicados a la codificación del régimen de giro 
    return genotipo_RPM
def obtener_binario_rc(individuo):
    genotipo_rc = individuo[8:12]       #Extrae los 4 bits dedicados a la codificación de la relación de compresión
    return genotipo_rc
def obtener_binario_masa_inyeccion(individuo):
    genotipo_masa = individuo[12:15]    #Extrae los 3 bits dedicados a la codificación de la masa de combustible secundario inyectado
    return genotipo_masa
def obtener_binario_composicion_combustible(individuo):
    #Extrae los 6 bits dedicados a la codificación de la presencia/ausencia de cada especie dentro del combustible inyectado
    presencia_componentes = individuo[15:22]    
    #Extrae los 12 bits dedicados a la codificación de la composición de cada especie dentro del combustible inyectado
    proporciones_componentes = individuo[22:]
    return presencia_componentes, proporciones_componentes

"""
Funciones para la descodificación de las cadenas de bits
"""
def decodificar_dosado(individuo):
    dosado_bits = individuo[:4]
    dosado_decimal = sum(bit * 2**i +0.001 for i, bit in enumerate(reversed(dosado_bits))) / 100 
    return dosado_decimal

def decodificar_rpm(individuo):
    rpm_bits = individuo[4:8]
    rpm_decimal = sum(bit * 2**i for i, bit in enumerate(reversed(rpm_bits)))
    rpm_decimal = int((rpm_decimal / 16) * 4000 + 1000)
    return rpm_decimal

def decodificar_relacion_compresion(individuo):
    rc_bits = individuo[8:12]
    rc_decimal = sum(bit * 2**i for i, bit in enumerate(reversed(rc_bits))) + 20
    return rc_decimal

def decodificar_masa_inyectada(individuo):
    masa_bits = individuo[12:15]
    masa_decimal = sum(bit * 2**i for i, bit in enumerate(reversed(masa_bits))) + 1
    return masa_decimal

def decodificar_composicion_inyeccion(individuo):
    ausencia_presencia = individuo[15:22]
    composicion = individuo[22:]
    # Almacena las composiciones decodificadas
    composicion_decodificada = []
    presencia_decodificada=dict(zip(especies_inyeccion,ausencia_presencia))

    # Decodifica la presencia/ausencia de cada especie y su composición
    for i in range(6):
        presencia = ausencia_presencia[i]
        composicion_bits = composicion[i*2:(i+1)*2]

        # Si la especie está presente, decodifica su composición
        if presencia == 1:
            if composicion_bits == [0, 0]:
                composicion_especie = 1
            elif composicion_bits == [0, 1]:
                composicion_especie = 2
            elif composicion_bits == [1, 0]:
                composicion_especie = 3
            else:
                composicion_especie = 4
        else:
            # Si la especie está ausente, asigna composición 0
            composicion_especie = 0

        # Agrega la composición decodificada a la lista
        composicion_decodificada.append(composicion_especie)
        
    suma=sum(composicion_decodificada)
    #Calcula la composición final atendiendo a las diferentes especies presentes
    composicion_final = {especie: valor / suma for especie, valor in zip(especies_inyeccion, composicion_decodificada)}

    return composicion_final




def check(Q,W,CO,CO2):
    """
    Función para comprobar la viabilidad de las soluciones.
    Toma como parámetros las variables cálculadas por la simulación
    """
    
    if W<=0:
        Apto = False
    elif W >= Q:
        Apto = False
    elif CO == 0 or CO2 == 0:
        Apto = False
    else: 
        Apto = True
    return Apto 


def Evaluacion_individuos(poblacion):
    
    """

    Función para la evaluación de los individuos presentes en la población de cada generación.
    
    Toma como parámetro poblacion (Lista que almacena los diferentes individuos a evaluar en cada generación)
    
    

    """
    matriz_puntuada = []
    valores_a_exportar=[]
    global historial_origen
    for i, individuo in  enumerate(poblacion):
        
        
        
        print(f"Índice: {i}")
        origen = next((o for (ind, o) in historial_origen if ind == individuo), None)
        print(f"Origen: {origen}")
        print("Individuo:", individuo)
        
        #Convierte el individuo a una tupla para comprobar si ya ha sido evaluado anteriormente
        #En caso afirmativo se recuperan las variables calculadas
        #En caso contrario se procede a evaluar el individuo
        individuo_tupla = tuple(individuo)
        if individuo_tupla in historico:
            print('Está en histórico')
            Fitness, _, Q, W, CO_emission, CO2_emission = historico[individuo_tupla]
            print('Fitness= ', Fitness)
            print(' ')
        else:
            
            #Variables a obtener del algoritmo
            dosado= decodificar_dosado(individuo)
            comp_inyeccion= decodificar_composicion_inyeccion(individuo)
            RPM= decodificar_rpm(individuo)                                          
            rc= decodificar_relacion_compresion(individuo)                               
            injector_mass =  decodificar_masa_inyectada(individuo)    
            
            
            #######################################################################################
            ###Definición de las variables del motor y de la simulación
            #######################################################################################
            
            
            inicio=tm.time()
            mecanismo='heptane-mehl.yaml'       #Mecanismo de reacción
            comp_air = 'o2:0.21, n2:0.79'       #Composición del aire
            comp_fuel = 'NC7H16:1'              #Composición del combustible principal
            air=ct.Solution('air.yaml')
            
            
            #Parámetros motor
            B=13e-2                                  #Diámetro
            S=16e-2                                  #Carrera
            c=rod_lenght=26.93e-2                    #
            a=S/2                                    #Longitud [m]
            A=0.25*np.pi*B**2                        #Área pistón
            
                                   
            f= RPM / 60                              #Frecuencia de giro [1/s]
            Vd = (np.pi*S* B**2)/ 4                  #Volumen desplazado [m^3]
            wdot=RPM*np.pi/30                        #Velocidad giro [rad/s]
            Vmax= Vd + Vd/(rc-1)                     #Volumen máximo
            
            #Parámetros Admisión 
            T_adm=500                                 #T admisión [ºK]
            P_adm= 1e5                                #P admisión [Pa]
            A_val_admision= -18/ 180. * np.pi         #Apertura válvula admisión [º]
            C_val_admision= 180/ 180. * np.pi         #Cierre válvula admisión [º]
            coef_valv_admision=1.e-6
            
            #Parámetros Escape
            T_escape= 350                             #T escape [ºK]
            P_escape= 1e5                             #P escape [Pa]
            A_val_escape= 522/ 180. * np.pi           #Apertura válvula admisión [º]
            C_val_escape= 20/ 180. * np.pi            #Cierre válvula admisión [º]
            coef_valv_escape=1.e-6
            
            # Inyección Combustible
            T_injector = 300                          #T inyección [ºK]
            p_injector = 1600e5                       #P inyección [Pa]
            comp_injector = comp_inyeccion            #Composición combustible a inyectar 
            injector_open = 350. / 180. * np.pi       #Apertura inyector [º]
            injector_close = 365. / 180. * np.pi      #Cierre inyección   [º]         
            injector_mass =  injector_mass*1e-5       #Masa combustible [Kg]
            
            # Parámetros Ambiente
            T_ambient = 300                           #T ambiente [ºK]
            p_ambient = 1e5                           #P ambiente [Pa]
            comp_ambient = comp_air                   #Composición Ambiente
            
            # Parámetros Simulación 
            sim_n_revolutions = 2
            delta_T_max = 20.
            rtol = 1.e-12
            atol = 1.e-16
            
            # Crea una lista para las especies a inyectar
            componentes = [f"{clave}:{valor}" for clave, valor in comp_inyeccion.items()]
            comp_inyeccion = ', '.join(componentes)
            comp_inyeccion= comp_inyeccion.replace(' ', '')
            
            
            
            #######################################################################################
            ###Definición Parámetros y Funciones motor HCCI 
            #######################################################################################
            
            def v_piston(t):
                """
                Cálcula la velocidad del pistón en función del instante de la simulación considerado
                """
                return  -(a**2*np.sin(wdot*t)*np.cos(wdot*t)*wdot)/np.sqrt(c**2-(a*np.sin(wdot*t))**2) -a*np.sin(wdot*t)*wdot 
            def Volume(t):
                """
                Cálcula el volumen en función del instante de la simulación considerado
                """
                return Vd/(rc-1)+ (0.25*np.pi*B**2)*(c+a-a*np.cos(wdot*t)-np.sqrt(c**2-(a*np.sin(wdot*t))**2))
            def angulo_giro(t):
                """
                Cálcula el ángulo de giro del cigüeñal en función del instante de la simulación considerado
                """
                return 180*wdot*t/np.pi
            def grados_a_rad(theta):
                """
                Convierte de grados a radianes
                """
                return theta*np.pi/180
            def obtener_composicion(comp_air, comp_fuel, dosado): 
                """
                Cálcula la composición en función de las composiciónes del aire y del combustible
                """
                comp_air = comp_air
                comp_fuel= comp_fuel
                dosado= dosado 
                fracciones_combustible = {}
                
                # Calcula las fracciones molares del aire
                frac_o2_air = float(comp_air.split(',')[0].split(':')[1])
                frac_n2_air = float(comp_air.split(',')[1].split(':')[1])
                
                for comp in comp_fuel.split(','): 
                    especie, fraccion = comp.split(':')
                    fracciones_combustible[especie]= float(fraccion)*dosado
                    
                # Calcula las fracciones molares ajustadas según la fracción deseada
                fraccion_air = 1 - dosado
                frac_o2 = frac_o2_air * fraccion_air
                frac_n2 = frac_n2_air * fraccion_air
                fracciones_ajustadas = {'O2': frac_o2, 'N2': frac_n2}
                for comp in fracciones_combustible:
                    fracciones_ajustadas[comp] = fracciones_combustible[comp]
              
                return fracciones_ajustadas
            
            #######################################################################################
            ###Definición Reactores
            #######################################################################################
            
            comp_admision= obtener_composicion(comp_air, comp_fuel, dosado)
            gas=ct.Solution(mecanismo)
            
            #Definición condiciones iniciales reactor
            gas.TPX=T_adm,P_adm,comp_admision
            r= ct.IdealGasReactor(gas, energy ='on')
            r.volume= Vd/(rc-1)
            
            
            ##### Válvula admisión  #####
            gas.TPX= T_adm,P_adm,comp_admision
            admision= ct.Reservoir(gas)
            admision_valve= ct.Valve(admision,r)
            admision_delta= np.mod(C_val_admision - A_val_admision, 4*np.pi )
            admision_valve.valve_coeff= coef_valv_admision
            
            # Definición de la función de tiempo para la válvula de admisión
            # Considera los ángulos de apertura y cierre de válvula
            admision_valve.time_function = (
               lambda t:np.mod( grados_a_rad(angulo_giro(t))- A_val_admision,  4 * np.pi) < admision_delta)
            
            ##### Inyector  #####
            gas.TPX = T_injector, p_injector, comp_injector
            injector = ct.Reservoir(gas)
            injector_mfc = ct.MassFlowController(injector, r)
            injector_delta = np.mod(injector_close - injector_open, 4 * np.pi)
            injector_t_open = (injector_close - injector_open) / 2. / np.pi / f
            injector_mfc.mass_flow_coeff = injector_mass / injector_t_open
            
            # Definición de la función de tiempo para el inyector
            # Considera los ángulos de apertura y cierre 
            injector_mfc.time_function = (
                lambda t: np.mod( grados_a_rad(angulo_giro(t)) - injector_open, 4 * np.pi) < injector_delta)
            
            ##### Válvula escape  #####
            gas.TP=T_escape, P_escape
            escape=ct.Reservoir(gas)
            escape_valve= ct.Valve(r, escape)
            escape_delta= np.mod( C_val_escape -  A_val_escape, 4*np.pi )
            escape_valve.valve_coeff= coef_valv_escape
            
            # Definición de la función de tiempo para la válvula de escape
            # Considera los ángulos de apertura y cierre de válvula
            escape_valve.time_function = (
                lambda t:np.mod( grados_a_rad(angulo_giro(t))- A_val_escape,  4 * np.pi) < escape_delta)
            
            ##### Ambiente  #####
            gas.TPX = T_ambient, p_ambient, comp_ambient
            ambient_air = ct.Reservoir(gas)
             
            ##### Pistón #####
            gas.TP = T_escape, P_escape
            env=ct.Reservoir(air)
            piston=ct.Wall(env,r,velocity=v_piston, A=A)
            
            ##### Simulación  #####
            sim=ct.ReactorNet([r])
            sim.rtol, sim.atol= rtol, atol
            
           #######################################################################################
           ### Ejecutar Simulación
           #######################################################################################
           
            states = ct.SolutionArray(
                r.thermo,
                extra=('t', 'ca', 'V', 'm', 'mdot_in', 'mdot_out', 'dWv_dt'),
                )
            
            dt = 1 / (360 * f)  
            t_stop = 1*sim_n_revolutions / f
            step_count=0
            while sim.time < t_stop:
                step_count+=1
                # Ajusta el paso de tiempo
                if angulo_giro(sim.time) > 180 and angulo_giro(sim.time) <= 360:
                    dt = 0.1 / (180 * f)       #Incrementa la resolución durante la carrera de compresión 
                    
                    
                # Realiza la integración en el tiempo 
                sim.advance(sim.time + dt)
               
                # Cálculo de las variables a guardar 
                dWv_dt = - (r.thermo.P - ambient_air.thermo.P) * A * \
                    v_piston(sim.time)
                
                #Almacena las variables calculadas durante la simulación
                states.append(r.thermo.state,
                              t=sim.time, ca=angulo_giro(sim.time),
                              V=r.volume, m=r.mass,
                              mdot_in=admision_valve.mass_flow_rate,
                              mdot_out=escape_valve.mass_flow_rate,
                              dWv_dt=dWv_dt)
            
            
            #Devuelve el tiempo empleado en la simulación
            fin = tm.time()
            print('Tiempo ejecución: ', -inicio+fin, 's')
            t = states.t
            
            
            ######################################################################
            ### Resultados integrales
            ######################################################################
            
            ##### Calor liberado #####
            Q = trapz(states.heat_release_rate * states.V, t)
            output_str = '{:45s}{:>4.2f} {}'
            print(output_str.format('Calor liberado:',
                                    Q / t[-1] / 1000., 'kW'))
            Q= Q / t[-1] / 1000
            
            ##### Trabajo #####
            W = trapz(states.dWv_dt, t)
            print(output_str.format('Trabajo de expansión:',
                                    W / t[-1] / 1000., 'kJ'))
            W= W / t[-1] / 1000.
            
            ##### Emisiones#####
            CO_emission = trapz(states('CO').mean_molecular_weight * states.mdot_out * states('CO').X[:, 0], t)
            CO_emission=CO_emission* 1.e6
            print(output_str.format('Emisiones CO:', CO_emission, 'ppm'))
            CO2_emission = trapz(states('CO2').mean_molecular_weight * states.mdot_out * states('CO2').X[:, 0], t)
            CO2_emission=CO2_emission* 1.e6
            print(output_str.format('Emisiones CO2:', CO2_emission, 'ppm'))
            
            
            
            ######################################################################
            ### Función de evaluación 
            ######################################################################
            
            if check(Q,W,CO_emission, CO2_emission):    #Comprueba que los individuos sean válidos
                
            
                #Coeficientes ponderación de la función de evaluación 
                A=2
                B=2
                C=8
                D=8
                epsilon=0.01
                
                #Normaliza las variables obtenidas
                Q_min=5.649845
                Q_max=69.317969
                W_min=0.181735
                W_max=18.511748
                CO_min=0.001236
                CO_max=5985.929725
                CO2_min=0.719296
                CO2_max=2124.459070
                
                F1_min=0.044030	
                F1_max=10.075916
                
                Q_norm=(Q-Q_min)/(Q_max-Q_min)
                W_norm=(W-W_min)/(W_max-W_min)
                CO_norm=(CO_emission-CO_min)/(CO_max-CO_min)
                CO2_norm=(CO2_emission-CO2_min)/(CO2_max-CO2_min)
                
                
                #Cálculo de las contribuciones
                F1= (W_norm + epsilon)/(Q_norm+ epsilon)
                F1_norm=(F1-F1_min)/(F1_max-F1_min)
                F2= W_norm
                F3= 0.01/ (CO_norm+ epsilon ) 
                F4=0.01/ (CO2_norm+ epsilon ) 
                
                Fitness=10*( A*F1_norm + B* F2 +C*F3 + D*F4 )
                
                print('Fitness=', Fitness)
                print('\n')
                
            else:                 #Asigna cero a todas las variables en caso de que el individuo no sea válido
                Fitness=0
                Q=0
                W=0
                CO_emission=0
                CO2_emission=0
                print('Individuo no válido')
                print('\n')
                
                #Se genera un nuevo individuo de la misma naturaleza que el individuo no válido
                
                origen = next((o for (i, o) in historial_origen if i == individuo), None) 
                if origen== 'Cruce':
                 nuevo_individuo= cruzamiento(seleccionados, 1, matriz_ordenada)  
                elif origen== 'Mutación':
                 nuevo_individuo= mutacion(seleccionados, tasa_mutacion, 1)
                else:
                 nuevo_individuo=creacion_individuos(1, 33, 6)
                 
                #Se añade el individuo creado a la población y al historial
                poblacion += nuevo_individuo
                historial_origen.append((nuevo_individuo[0], origen))
        
        #Se añade el individuo evaluado y las variables a matriz_puntuada
        matriz_puntuada.append((individuo, Fitness,Q,W,CO_emission,CO2_emission))
    
    return matriz_puntuada
    

def ordenar_matriz_puntuacion(matriz_puntuada):
    """
    Función que ordena la matriz_puntuada según la puntuación de los individuos
    """
    matriz_ordenada= sorted(matriz_puntuada, key= lambda x: x[1], reverse= True)
    return matriz_ordenada


#####################################################################################
### Funciones para la selección de los individuos 
#####################################################################################


def seleccion(matriz_ordenada, num_seleccionados):
    """
    Función para la selección de individuos a reproducirse

    Método de Selección: ELITISTA
    """
    seleccionados = []
    individuos_unicos = set()

    for tupla in matriz_ordenada:
        individuo = tupla[0]
        fitness = tupla[1]
        individuo_tupla = tuple(individuo)  # Convierte la lista en una tupla para que pueda ser añadida a un set
        if individuo_tupla not in individuos_unicos:
            seleccionados.append((individuo))
            individuos_unicos.add(individuo_tupla)
        if len(seleccionados) == num_seleccionados:
            break

    return seleccionados

def seleccion_por_ruleta(matriz_puntuada, num_seleccionados):
    """
    Función para la selección de individuos a reproducirse

    Método de Selección: RULETA
    """
    # Calcula la suma total de los puntajes de la población
    suma_total = sum(puntaje for _, puntaje in matriz_puntuada)
    
    # Genera una lista de probabilidades acumuladas 
    probabilidades_acumuladas = [sum(puntaje for _, puntaje in matriz_puntuada[:i+1]) / suma_total for i in range(len(matriz_puntuada))]
    
    # Selecciona individuos utilizando la ruleta
    seleccionados = []
    for _ in range(num_seleccionados):
        r = random.random()  # Generar un número aleatorio entre 0 y 1
        # Encontrar el individuo cuya probabilidad acumulada sea la más cercana a r
        for i, probabilidad_acumulada in enumerate(probabilidades_acumuladas):
            if r <= probabilidad_acumulada:
                seleccionados.append(matriz_puntuada[i][0])  
                break
    
    return seleccionados

def seleccion_por_torneo(matriz_puntuada, num_seleccionados, tam_torneo):
    """
    Función para la selección de individuos a reproducirse

    Método de Selección: TORNEO
    """
    seleccionados = []
    for _ in range(num_seleccionados):
        # Selecciona aleatoriamente un subconjunto de individuos para el torneo
        participantes = random.sample(matriz_puntuada, tam_torneo)
        # Selecciona al mejor individuo dentro del torneo
        mejor_individuo = max(participantes, key=lambda x: x[1])
        seleccionados.append(mejor_individuo[0])  # Agregar el individuo seleccionado
        
    return seleccionados


    
def partir_cromosomas(individuo):
    """
    Segmenta la cadena de bits del individuo en los diferentes genes
    """
    cromosoma1 = individuo[:4]
    cromosoma2 = individuo[4:8]
    cromosoma3 = individuo[8:12]
    cromosoma4 = individuo[12:15]
    cromosoma5 = individuo[15:]
    return cromosoma1, cromosoma2, cromosoma3, cromosoma4, cromosoma5

todos_los_hijos = []

def cruzamiento(seleccionados, num_descendencia, matriz_ordenada):
    """
    Función para realizar el cruzamiento entre los individuos seleccionados de la población
    Toma como parámetros la lista de individuos seleccionados, el número de hijos deseados y la matriz ordenada.
    """
    
    global todos_los_hijos
    hijos_cruce = []
    candidatos = seleccionados.copy()
    matriz_idx = len(seleccionados)

    while num_descendencia > 0:
        if len(candidatos) < 2:  #Si hay menos de dos padres se para el bucle
            break
        
        #Se elige al azar dos padres del conjunto de seleccionados
        padre1 = random.choice(candidatos)
        padre2 = random.choice(candidatos)
        
        # Asegurarse de que padre1 y padre2 sean diferentes
        while padre1 == padre2:
            padre2 = random.choice(candidatos)
        
        # Realizar cruzamiento. Llama a la función realizar_cruzamiento2
        hijo = realizar_cruzamiento2(padre1, padre2, matriz_ordenada, matriz_idx)
        
        #Si el hijo obtenido es nuevo se añade a los hijos creados, así como al histórico de hijos creados
        if hijo not in todos_los_hijos:
            hijos_cruce.append(hijo)
            todos_los_hijos.append(hijo)
            num_descendencia -= 1  # Decrementa el contador de descendencia

    return hijos_cruce

def realizar_cruzamiento2(padre1, padre2, matriz_ordenada, matriz_idx, max_intentos=100):
    """
    Función para llevar a cabo el cruzamiento entre los padres seleccionados. 
    Se realiza el cruzamiento mediante la asignación aleatoria de los genes de los padres al individuo nuevo.
    """
    
    global todos_los_hijos, historico
    cromosoma11, cromosoma21, cromosoma31, cromosoma41, cromosoma51 = partir_cromosomas(padre1)
    cromosoma12, cromosoma22, cromosoma32, cromosoma42, cromosoma52 = partir_cromosomas(padre2)

    intentos = 0
    while intentos < max_intentos:
        
        #Se asignan al azar los genes de los padres 
        cromosoma13 = random.choice([cromosoma11, cromosoma12])
        cromosoma23 = random.choice([cromosoma21, cromosoma22])
        cromosoma33 = random.choice([cromosoma31, cromosoma32])
        cromosoma43 = random.choice([cromosoma41, cromosoma42])
        cromosoma53 = random.choice([cromosoma51, cromosoma52])
        hijo = cromosoma13 + cromosoma23 + cromosoma33 + cromosoma43 + cromosoma53
        hijo_tupla = tuple(hijo)

        # Verifica si el hijo es diferente a los padres y no está ni en la lista global de hijos ni en el histórico
        if hijo != padre1 and hijo != padre2 and hijo not in todos_los_hijos and hijo_tupla not in historico:
            return hijo
        
        intentos += 1

    # Si no se puede generar un nuevo hijo después del número máximo de intentos, se usa el siguiente mejor individuo
    while matriz_idx < len(matriz_ordenada):
        siguiente_individuo = matriz_ordenada[matriz_idx][0]
        matriz_idx += 1
        
        #La función se llama a sí misma para realizar este nuevo cruzamiento
        hijo = realizar_cruzamiento2(padre1, siguiente_individuo, matriz_ordenada, matriz_idx, max_intentos)
        if hijo not in todos_los_hijos:
            return hijo

    # Si no se puede generar un nuevo hijo devuelve uno de los padres
    return random.choice([padre1, padre2])


def comprobar_individuos(individuo):
    """
    Función para comprobar que los individuos mutados cumplen con los requisitos impuestos sobre la presencia de especies
    """
    P1 = individuo[15]
    P2 = individuo[16]
    P3 = individuo[17]
    P4 = individuo[18]
    P5 = individuo[19]
    P6 = individuo[20]
    suma= P1 + P2 + P3+ P4+ P5 +P6
    if  suma < 1 or suma > 3:
        return False
    if P1 == 0 and individuo[21:23] != [0, 0]: 
        return False
    if P2 == 0 and individuo[23:25] != [0, 0]: 
        return False
    if P3 == 0 and individuo[25:27] != [0, 0]: 
        return False
    if P4 == 0 and individuo[27:29] != [0, 0]: 
        return False
    if P5 == 0 and individuo[29:31] != [0, 0]: 
        return False
    if P6 == 0 and individuo[31:33] != [0, 0]: 
        return False
        
    return True


todos_los_individuos_mutados = []

def mutacion(individuos, tasa_mutacion, n_a_mutar):
    """
    Función para realizar la mutación de los individuos.
    Recibe como parámetros los posibles individuos a mutar, la tasa de mutación y el número de individuos que se desea mutar
    La mutación se realiza modificando los bits: 0 a 1 o 1 a 0
    """
    global todos_los_individuos_mutados
    individuos_mutados = []
    individuos_a_mutar = random.sample(individuos, n_a_mutar)  # Selecciona al azar n_a_mutar individuos
                                                                
    for individuo in individuos_a_mutar:
        intentos = 0
        max_intentos = 200  
        while intentos < max_intentos:
            individuo_mutado = list(individuo)  # Convierte al individuo en una lista mutable
            for i in range(len(individuo_mutado) - 18):
                if random.random() < tasa_mutacion:  # Comprueba si se debe mutar este bit
                    individuo_mutado[i] = 1 - individuo_mutado[i]  # Mutación: cambiar 0 a 1 o 1 a 0
            
            if comprobar_individuos(individuo_mutado) and individuo_mutado not in todos_los_individuos_mutados:
                individuos_mutados.append(individuo_mutado)
                todos_los_individuos_mutados.append(individuo_mutado)
                break
            intentos += 1
        
        if intentos == max_intentos:  #Si se alcanza el número máximo de intentos se selecciona otro individuo a mutar 
            nuevo_individuo = random.choice(individuos)
            # Realiza la mutación en el nuevo individuo seleccionado
            nuevo_individuo_mutado = mutacion([nuevo_individuo], tasa_mutacion, 1)[0]
            individuos_mutados.append(nuevo_individuo_mutado)
            todos_los_individuos_mutados.append(nuevo_individuo_mutado)
    
    return individuos_mutados

def modificar_poblacion(poblacion, seleccionados, hijos_cruce, individuos_mutados, nuevo_individuo):
    """
    Función para realizar la actualización de la población
    Recibe como parámetros la población actual, los individuos seleccionados, los individuos obtenidos del cruzamiento, 
    los individuos procedentes de la mutación y finalmente los nuevos individuos creados.
    """
    # Se añaden los individuos seleccionados
    nueva_poblacion = list(seleccionados)
    
    # Se añaden los individuos resultantes del cruce
    nueva_poblacion.extend(hijos_cruce)
    
    # Se añaden los individuos mutados
    nueva_poblacion.extend(individuos_mutados)
    
    # Se añaden los individuos creados para completar la población
    nueva_poblacion.extend(nuevo_individuo)
    
    #Se actualiza el historial de origen para realizar el seguimiento de los individuos
    for individuo in seleccionados:
        historial_origen.append((individuo,'Selección'))
        
    for _ in hijos_cruce:
        historial_origen.append((_,'Cruce'))
   
    for _ in individuos_mutados:
        historial_origen.append((_,'Mutación'))
    
    for _ in nuevo_individuo:
        historial_origen.append((_,'Creación'))
       
    return nueva_poblacion




#####################################################################################
### Desarrollo del Algoritmo Genético
#####################################################################################

    
maximo_generaciones=10          # Máximo de generaciones a estudiar

numero_de_generaciones =0       # Inicializa el contador de generaciones
num_individuos =10              # Número de individuos a crear en la primera generación


num_seleccionados=2             # Número de individuos a seleccionar
num_descendencia=1              # Número de hijos a obtener
n_a_mutar0=2                    # Número de individuos a mutar
tasa_mutacion=0.15              # Tasa de mutación
num_individuos_a_crear=2        # Número de individuos a crear aleatoriamente



#### Inicialización de variables ####

poblacion = []
historico={}
historico_poblacion=[]   

#### Inicialización de la población  ####

individuos = creacion_individuos(num_individuos, 33, 6)                         #Se generan los individuos
historial_origen = [(individuo, 'Creación') for individuo in individuos]        #Se actualiza el historial de origen
poblacion += individuos



#####################################################################################
### Bucle Principal
#####################################################################################

while numero_de_generaciones < maximo_generaciones:
    n_a_mutar=n_a_mutar0
    print("Generación:", numero_de_generaciones)
    
    # Evaluación individuos 
    matriz_puntuada = Evaluacion_individuos(poblacion)
    
    # Se ordenan los individuos por puntuación                      
    matriz_ordenada = ordenar_matriz_puntuacion(matriz_puntuada)       
       
    # Se añaden al histórico los individuos nuevos
    for tupla in matriz_ordenada:                                             
        individuo = tupla[0]
        puntuacion = tupla[1]
        Q=tupla[2]
        W=tupla[3]
        CO_emission=tupla[4]
        CO2_emission=tupla[5]
        individuo_tupla = tuple(individuo)  # Convertimos la lista en una tupla, para poder acceder a los términos
        if individuo_tupla not in historico: 
            historico[individuo_tupla] = (puntuacion, numero_de_generaciones ,Q, W, CO_emission, CO2_emission)
    historico_poblacion.append((matriz_puntuada)) 
    
    # Selección individuos
    seleccionados = seleccion(matriz_ordenada, num_seleccionados) 
    
    #Cruzamiento
    hijos_cruce = cruzamiento(seleccionados, num_descendencia, matriz_ordenada)
    if len(hijos_cruce) < num_descendencia:         #Asegura que en caso de que no haya sido posible obtener la descendencia deseada mediante cruzamiento se genere un nuevo individuo, en este caso, por mutación
        n_a_mutar+= 1
        
    #Mutación
    mutados = mutacion(seleccionados, tasa_mutacion, n_a_mutar)
    
    #Creación nuevos individuos
    nuevo_individuo = creacion_individuos(num_individuos_a_crear, 33, 6)
    
    #Actualización Población
    poblacion = modificar_poblacion(poblacion, seleccionados, hijos_cruce, mutados, nuevo_individuo)

    # Imprimir la mejor puntuación y la media de la puntuación de la generación actual
    mejor_puntuacion = matriz_ordenada[0][1]
    print("Mejor puntuación de la generación:", mejor_puntuacion)
    matriz_limpia = [tupla for tupla in matriz_ordenada if tupla[1] != 0.0]
    puntuaciones  = [tupla[1] for tupla in matriz_limpia]
    media_puntuaciones = sum(puntuaciones) / len(puntuaciones)
    print("Media de la puntuación de la generación:", media_puntuaciones)
    print('Mejor Individuo:', matriz_ordenada[0][0])
    print(' ')

    # Incrementa el contador de generaciones
    numero_de_generaciones += 1
    
    # Para evitar fallos al guardar los archivos pongo la fecha con hora en el nombre 
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")  # Formato: YYYYMMDD_HHMMSS

    # Verifica la condición de parada 
    if numero_de_generaciones >= maximo_generaciones:
        
        #####################################################################################
        ### EXPORTACIONES
        #####################################################################################

        #Nombre columnas a exportar
        columnas = ["Individuos", "Fitness", "Generación", "Calor Liberado", "Trabajo de Expansión", "Emisiones de CO", "Emisiones de CO2"]
    
        
        now = datetime.datetime.now()
        timestamp = now.strftime("%Y%m%d_%H%M%S")  
        
        # Nombres de los archivos CSV a exportar
        nombre_prueba_csv = f"Historicocaso4_{timestamp}.csv"
        nombre_origen_csv = f"Origencaso4_{timestamp}.csv"
        
        # Exportar datos del historico
        with open(nombre_prueba_csv, "w", newline='') as archivo_csv:
            escritor_csv = csv.writer(archivo_csv)
            archivo_csv.write(f"# Máximo de generaciones: {maximo_generaciones}, Número de generaciones: {numero_de_generaciones}, Número de individuos: {num_individuos},\n")
            archivo_csv.write(f"# Número de seleccionados: {num_seleccionados}, Número de descendencia: {num_descendencia}, Número a mutar: {n_a_mutar}, Tasa de mutación: {tasa_mutacion}, Número de individuos a crear: {num_individuos_a_crear}\n")
    
            
            # Escribe los encabezados de las columnas
            escritor_csv.writerow(columnas)
            
            # Escribe los datos de la matriz ordenada
            for generacion_idx, generacion in enumerate(historico_poblacion):
                for individuo in generacion:
                    cromosoma, fitness, Q, W, CO_emission, CO2_emission = individuo
                    fila = [cromosoma, fitness, generacion_idx, Q, W, CO_emission, CO2_emission]
                    escritor_csv.writerow(fila)
                
        # Exportar datos del historial_origen 
        with open(nombre_origen_csv, "w", newline='') as archivo_csv:
            escritor_csv = csv.writer(archivo_csv)
            
            escritor_csv.writerow(["Historial Origen"])
            
            # Escribe los datos de historial_origen 
            for valor in historial_origen:
                escritor_csv.writerow([valor])
