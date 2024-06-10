# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 19:55:39 2024

@author: Pablo Ramos
"""

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


def creacion_individuos(num_individuos, longitud, especies):
    """
    Función para la creación de nuevos individuos. 
    -num_individuos indica el tamaño de la población a crear 
    -longitud indica el tamaño de la cadena binaria para cada uno de los individuos
    
    La variable longitud vendrá determinada principalmente por las especies contempladas 
    """
    individuos_unicos = set()
    nuevos_individuos = []
    
    # Mientras no tengamos suficientes individuos únicos, seguimos generando nuevos individuos
    while len(individuos_unicos) < num_individuos:
        nuevo_individuo = tuple(np.random.randint(2, size=longitud-2*especies))  # Generar un nuevo individuo como una tupla de bits
        presencia_especies = nuevo_individuo[15:] 
        
        # Calcular la suma de los bits de presencia/ausencia
        suma = sum(presencia_especies)
        
        # Verificar que el número de especies presentes sea menor o igual a 3 
        if 1 <= suma <= 3:
            #Asignamos los bits de composición según los bits de presencia/ausencia
            for i, presencia in enumerate(presencia_especies):
                if presencia == 1: 
                    composicion = tuple(np.random.randint(2, size=2))
                    nuevo_individuo += composicion
                else:
                    composicion = (0, 0)
                    nuevo_individuo += composicion
                    
            individuos_unicos.add(nuevo_individuo)
            # Agregar el nuevo individuo al conjunto de individuos únicos
    
    # Convertir el conjunto de individuos únicos de nuevo a una lista de listas de bits
    for individuo in individuos_unicos:
        nuevos_individuos.append(list(individuo))
        
    return nuevos_individuos


"""
Funciones para segmentar la cadena de bits del individuo en función de la variable a la que representan 
"""
def obtener_binario_dosado(individuo):
    genotipo_dosado = individuo[0:4]
    return genotipo_dosado
def obtener_binario_RPM(individuo):
    genotipo_RPM = individuo[4:8]
    return genotipo_RPM
def obtener_binario_rc(individuo):
    genotipo_rc = individuo[8:12]
    return genotipo_rc
def obtener_binario_masa_inyeccion(individuo):
    genotipo_masa = individuo[12:15]
    return genotipo_masa
def obtener_binario_composicion_combustible(individuo):
    # Supongamos que los primeros 6 bits representan la presencia o ausencia de cada componente
    presencia_componentes = individuo[15:22]
    # Supongamos que los siguientes 12 bits representan las proporciones de cada componente
    proporciones_componentes = individuo[22:]
    return presencia_componentes, proporciones_componentes

"""
Funciones para la descodificación de las cadenas de bits
"""
def decodificar_dosado(individuo):
    dosado_bits = individuo[:4]
    dosado_decimal = sum(bit * 2**i +0.001 for i, bit in enumerate(reversed(dosado_bits))) / 100 #15.0
    #print('Dosado:', dosado_decimal)
    return dosado_decimal

def decodificar_rpm(individuo):
    rpm_bits = individuo[4:8]
    rpm_decimal = sum(bit * 2**i for i, bit in enumerate(reversed(rpm_bits)))
    rpm_decimal = int((rpm_decimal / 16) * 4000 + 1000)
    #print('RPM:', rpm_decimal)
    return rpm_decimal

def decodificar_relacion_compresion(individuo):
    rc_bits = individuo[8:12]
    rc_decimal = sum(bit * 2**i for i, bit in enumerate(reversed(rc_bits))) + 20
    #print('RC:', rc_decimal)
    return rc_decimal

def decodificar_masa_inyectada(individuo):
    masa_bits = individuo[12:15]
    masa_decimal = sum(bit * 2**i for i, bit in enumerate(reversed(masa_bits))) + 1
    #print('Masa inyección:',masa_decimal)
    return masa_decimal

def decodificar_composicion_inyeccion(individuo):
    ausencia_presencia = individuo[15:22]
    composicion = individuo[22:]
    #print(ausencia_presencia)
    #print(composicion)

    # Lista para almacenar las composiciones decodificadas
    composicion_decodificada = []
    presencia_decodificada=dict(zip(especies_inyeccion,ausencia_presencia))

    # Decodificar la presencia/ausencia de cada especie y su composición
    for i in range(6):
        presencia = ausencia_presencia[i]
        composicion_bits = composicion[i*2:(i+1)*2]

        # Si la especie está presente, decodificar su composición
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
            # Si la especie está ausente, asignar composición 0
            composicion_especie = 0

        # Agregar la composición decodificada a la lista
        composicion_decodificada.append(composicion_especie)
        
    suma=sum(composicion_decodificada)
    composicion_final = {especie: valor / suma for especie, valor in zip(especies_inyeccion, composicion_decodificada)}

    return composicion_final




"""
Función para comprobar la viabilidad de las soluciones
"""

def check(Q,W,CO,CO2):
    
    if W<=0:
        Apto = False
    elif W >= Q:
        Apto = False
    elif CO == 0 or CO2 == 0:
        Apto = False
    else: 
        Apto = True
    return Apto 

"""
Función para obtención de los gráficos
"""
def obtener_graficas(states, comp_fuel, comp_inyeccion):
#Gráfico de concentraciones de especies
        plt.figure(figsize=(10, 6))
        
        especies_graficar=['O2','CO','CO2']
        for comp in comp_fuel.split(','): 
            especie = comp.split(':')[0]
            especies_graficar.append(especie)
        for comp in comp_inyeccion.split(','): 
            especie = comp.split(':')[0]
            especies_graficar.append(especie)
        for especie in especies_graficar:
            plt.plot(states.ca, states(especie).X, label=especie)
        
        plt.legend(loc='upper right')
        plt.ylabel('$X_i$ [-]')
        plt.xlabel(r'$\phi$ [º]')
        plt.gca().xaxis.set_major_locator(MultipleLocator(180))
        plt.title('Evolución de Concentraciones de Especies en función del ángulo del cigüeñal')
        plt.grid(True)
        plt.show()
        
        # Gráfico de volumen
        #plt.figure(figsize=(14, 10))
        #plt.plot(states.ca, states.V)
        #plt.ylabel('Volumen')
        #plt.gca().xaxis.set_major_locator(MultipleLocator(180))
        #plt.xlabel('Ángulo del Cigüeñal (grados)')
        #plt.title('Evolución del Volumen en función del ángulo del cigüeñal')
        #plt.grid(True)
        #plt.show()
        
        # Gráfico de presión
        plt.figure(figsize=(10, 6))
        plt.plot(states.ca, states.P)
        plt.ylabel('Presión (Pa)')
        plt.xlabel('Ángulo del Cigüeñal (grados)')
        plt.gca().xaxis.set_major_locator(MultipleLocator(180))
        plt.title('Evolución de la Presión en función del ángulo del cigüeñal')
        plt.grid(True)
        plt.show()
        
        # Gráfico de temperatura
        plt.figure(figsize=(10, 6))
        plt.plot(states.ca, states.T)
        plt.ylabel('Temperatura (K)')
        plt.xlabel('Ángulo del Cigüeñal (grados)')
        plt.gca().xaxis.set_major_locator(MultipleLocator(180))
        plt.title('Evolución de la Temperatura en función del ángulo del cigüeñal')
        plt.grid(True)
        plt.show()
        
        # Gráfico de flujo másico por la válvula de escape y entrada
        plt.figure(figsize=(10, 6))
        plt.plot(states.ca, states.mdot_in, label='Flujo másico de entrada')
        plt.plot(states.ca, states.mdot_out, label='Flujo másico de salida')
        plt.legend(loc='upper right')
        plt.ylabel('Flujo másico (kg/s)')
        plt.xlabel('Ángulo del Cigüeñal (grados)')
        plt.gca().xaxis.set_major_locator(MultipleLocator(180))
        plt.title('Flujo másico a través de las válvulas en función del ángulo del cigüeñal')
        plt.grid(True)
        plt.show()
    
def Evaluacion_individuos(poblacion):
    matriz_puntuada = []
    valores_a_exportar=[]
    for i, individuo in  enumerate(poblacion):
        
        
        
        print(f"Índice: {i}")
        print("Individuo:", individuo)
        
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
            RPM= decodificar_rpm(individuo)                                             #Revoluciones por minuto
            rc= decodificar_relacion_compresion(individuo)                              #Relación de compresión 
            injector_mass =  decodificar_masa_inyectada(individuo)    
            
            
            #Combustible y aire 
            inicio=tm.time()
            mecanismo='heptane-mehl.yaml'
            comp_air = 'o2:0.21, n2:0.79'
            comp_fuel = 'NC7H16:1'          # f'ch4:{0.8*dosado},c2h6:{0.1*dosado},c3h8:{0.1*dosado}'
            air=ct.Solution('air.yaml')
            
            
            #Parámetros motor
            B=13e-2
            S=16e-2
            c=rod_lenght=26.93e-2
            a=S/2                               #Longitud [m]
            A=0.25*np.pi*B**2
            
                                   
            f= RPM / 60                         #Régimen de giro [1/s]
            Vd = (np.pi*S* B**2)/ 4             #Volumen desplazado [m^3]
            wdot=RPM*np.pi/30                   #Velocidad giro [rad/s]
            Vmax= Vd + Vd/(rc-1)         
            
            #Parámetros Admisión 
            T_adm=500
            P_adm= 1e5
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
            comp_injector = comp_inyeccion
            injector_open = 350. / 180. * np.pi       #Apertura inyector [º]
            injector_close = 365. / 180. * np.pi      #Cierre inyección   [º]         
            injector_mass =  injector_mass*1e-5       #Masa combustible [Kg]
            
            # Parámetros Ambiente
            T_ambient = 300.                          #T ambiente [ºK]
            p_ambient = 1e5                           #P ambiente [Pa]
            comp_ambient = comp_air
            
            # Parámetros Simulación 
            sim_n_revolutions = 2
            delta_T_max = 20.
            rtol = 1.e-12
            atol = 1.e-16
            
            # Crear una lista para las especies a inyectar
            componentes = [f"{clave}:{valor}" for clave, valor in comp_inyeccion.items()]
            comp_inyeccion = ', '.join(componentes)
            comp_inyeccion= comp_inyeccion.replace(' ', '')
            
            
            #####Definición Parámetros y Funciones motor HCCI ############################
            
            
            def v_piston(t):
                return  -(a**2*np.sin(wdot*t)*np.cos(wdot*t)*wdot)/np.sqrt(c**2-(a*np.sin(wdot*t))**2) -a*np.sin(wdot*t)*wdot 
            def Volume(t):
                return Vd/(rc-1)+ (0.25*np.pi*B**2)*(c+a-a*np.cos(wdot*t)-np.sqrt(c**2-(a*np.sin(wdot*t))**2))
            def angulo_giro(t): 
                return 180*wdot*t/np.pi
            def grados_a_rad(theta):
                return theta*np.pi/180
            def obtener_composicion(comp_air, comp_fuel, dosado): 
                
                comp_air = comp_air
                comp_fuel= comp_fuel
                dosado= dosado 
                fracciones_combustible = {}
                
                # Calcular las fracciones molares del aire
                frac_o2_air = float(comp_air.split(',')[0].split(':')[1])
                frac_n2_air = float(comp_air.split(',')[1].split(':')[1])
                
                for comp in comp_fuel.split(','): 
                    especie, fraccion = comp.split(':')
                    fracciones_combustible[especie]= float(fraccion)*dosado
                    
                # Calcular las fracciones molares ajustadas según la fracción deseada
                fraccion_air = 1 - dosado
                frac_o2 = frac_o2_air * fraccion_air
                frac_n2 = frac_n2_air * fraccion_air
                fracciones_ajustadas = {'O2': frac_o2, 'N2': frac_n2}
                for comp in fracciones_combustible:
                    fracciones_ajustadas[comp] = fracciones_combustible[comp]
              
                return fracciones_ajustadas
            
            ######Definición Reactores####################################################
            
            comp_admision= obtener_composicion(comp_air, comp_fuel, dosado)
            
            gas=ct.Solution(mecanismo)
            #Definición condiciones iniciales reactor
            gas.TPX=T_adm,P_adm,comp_admision
            r= ct.IdealGasReactor(gas, energy ='on')
            r.volume= Vd/(rc-1)
            
            
            #Válvula admisión 
            gas.TPX= T_adm,P_adm,comp_admision
            admision= ct.Reservoir(gas)
            admision_valve= ct.Valve(admision,r)
            admision_delta= np.mod(C_val_admision - A_val_admision, 4*np.pi )
            admision_valve.valve_coeff= coef_valv_admision
            # Definición de la función de tiempo para la válvula de admisión
            admision_valve.time_function = (
               lambda t:np.mod( grados_a_rad(angulo_giro(t))- A_val_admision,  4 * np.pi) < admision_delta)
            
            #Inyector
            gas.TPX = T_injector, p_injector, comp_injector
            injector = ct.Reservoir(gas)
            
            injector_mfc = ct.MassFlowController(injector, r)
            injector_delta = np.mod(injector_close - injector_open, 4 * np.pi)
            injector_t_open = (injector_close - injector_open) / 2. / np.pi / f
            injector_mfc.mass_flow_coeff = injector_mass / injector_t_open
            injector_mfc.time_function = (
                lambda t: np.mod( grados_a_rad(angulo_giro(t)) - injector_open, 4 * np.pi) < injector_delta)
            
            #Válvula escape 
            gas.TP=T_escape, P_escape
            escape=ct.Reservoir(gas)
            escape_valve= ct.Valve(r, escape)
            escape_delta= np.mod( C_val_escape -  A_val_escape, 4*np.pi )
            escape_valve.valve_coeff= coef_valv_escape
            escape_valve.time_function = (
                lambda t:np.mod( grados_a_rad(angulo_giro(t))- A_val_escape,  4 * np.pi) < escape_delta)
            
            #Ambiente
            gas.TPX = T_ambient, p_ambient, comp_ambient
            ambient_air = ct.Reservoir(gas)
             
            #Pistón 
            gas.TP = T_escape, P_escape
            env=ct.Reservoir(air)
            piston=ct.Wall(env,r,velocity=v_piston, A=A)
            
            #Simulación
            sim=ct.ReactorNet([r])
            sim.rtol, sim.atol= rtol, atol
            
            
            ###Ejecutar Simulación#######################################################
            
            states = ct.SolutionArray(
                r.thermo,
                extra=('t', 'ca', 'V', 'm', 'mdot_in', 'mdot_out', 'dWv_dt'),
                )
            
            dt = 1 / (360 * f)  
            t_stop = 1*sim_n_revolutions / f
            step_count=0
            while sim.time < t_stop:
                #print(f"Time: {sim.time}, Steps: {step_count}")
                step_count+=1
                # Ajusta el paso de tiempo durante la carrera de compresión
                if angulo_giro(sim.time) > 180 and angulo_giro(sim.time) <= 360:
                    dt = 0.1 / (180 * f)                      #Incrementa la resolución durante la carrera de compresión 
                    #print('Time: {}, T: {:.2f} ºK, P:{:.3f} bar'.format(sim.time,r.thermo.T, r.thermo.P*1e-5))
                    
                # Realiza la integración en el tiempo 
                sim.advance(sim.time + dt)
               
                # Cálculo de las variables a guardar 
                dWv_dt = - (r.thermo.P - ambient_air.thermo.P) * A * \
                    v_piston(sim.time)
                    
                states.append(r.thermo.state,
                              t=sim.time, ca=angulo_giro(sim.time),
                              V=r.volume, m=r.mass,
                              mdot_in=admision_valve.mass_flow_rate,
                              mdot_out=escape_valve.mass_flow_rate,
                              dWv_dt=dWv_dt)
            
            fin = tm.time()
            print('Tiempo ejecución: ', -inicio+fin, 's')
            #print('\n')
            t = states.t
            
            
            ######################################################################
            # Resultados integrales
            ######################################################################
            
            # Calor liberado
            Q = trapz(states.heat_release_rate * states.V, t)
            output_str = '{:45s}{:>4.2f} {}'
            print(output_str.format('Calor liberado:',
                                    Q / t[-1] / 1000., 'kW'))
            Q= Q / t[-1] / 1000
            
            # Trabajo 
            W = trapz(states.dWv_dt, t)
            print(output_str.format('Trabajo de expansión:',
                                    W / t[-1] / 1000., 'kJ'))
            W= W / t[-1] / 1000.
            # Emisiones 
            CO_emission = trapz(states('CO').mean_molecular_weight * states.mdot_out * states('CO').X[:, 0], t)
            CO_emission=CO_emission* 1.e6
            print(output_str.format('Emisiones CO:', CO_emission, 'ppm'))
            CO2_emission = trapz(states('CO2').mean_molecular_weight * states.mdot_out * states('CO2').X[:, 0], t)
            CO2_emission=CO2_emission* 1.e6
            print(output_str.format('Emisiones CO2:', CO2_emission, 'ppm'))
            
            
            
            ######################################################################
            # Función de evaluación 
            ######################################################################
            if check(Q,W,CO_emission, CO2_emission):
                
                A=7
                B=5.5
                C=0.5
                D=0.5
                epsilon=0.01
                
              
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
                
                F1= (W_norm + epsilon)/(Q_norm+ epsilon)
                F1_norm=(F1-F1_min)/(F1_max-F1_min)
                F2= W_norm
                F3= 0.01/ (CO_norm+ epsilon ) 
                F4=0.01/ (CO2_norm+ epsilon ) 
                
                Fitness=10*( A*F1_norm + B* F2 +C*F3 + D*F4 )
                
                #matriz_anexa[fila,0]= Fitness
                #fila+= 1
                print('Fitness=', Fitness)
                print('\n')
            else: 
                Fitness=0
                Q=0
                W=0
                CO_emission=0
                CO2_emission=0
                print('Individuo no válido')
                print('\n')
                nuevo_individuo=creacion_individuos(1, 33, 6)
                historial_origen.append('Creación')
                poblacion += nuevo_individuo
            
        matriz_puntuada.append((individuo, Fitness,Q,W,CO_emission,CO2_emission))
    
    return matriz_puntuada
    
def ordenar_matriz_puntuacion(matriz_puntuada):
    matriz_ordenada= sorted(matriz_puntuada, key= lambda x: x[1], reverse= True)
    return matriz_ordenada


#####################################################################################
# Funciones para la selección de los individuos 
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
        individuo_tupla = tuple(individuo)  # Convertir la lista en una tupla para que pueda ser añadida a un set
        if individuo_tupla not in individuos_unicos:
            seleccionados.append((individuo))
            individuos_unicos.add(individuo_tupla)
        if len(seleccionados) == num_seleccionados:
            break

    return seleccionados

def seleccion_por_ruleta(matriz_puntuada, num_seleccionados):
    # Calcular la suma total de los puntajes de la población
    suma_total = sum(puntaje for _, puntaje in matriz_puntuada)
    
    # Generar una lista de probabilidades acumuladas para la selección
    probabilidades_acumuladas = [sum(puntaje for _, puntaje in matriz_puntuada[:i+1]) / suma_total for i in range(len(matriz_puntuada))]
    
    # Seleccionar individuos utilizando la ruleta
    seleccionados = []
    for _ in range(num_seleccionados):
        r = random.random()  # Generar un número aleatorio entre 0 y 1
        # Encontrar el individuo cuya probabilidad acumulada sea la más cercana a r
        for i, probabilidad_acumulada in enumerate(probabilidades_acumuladas):
            if r <= probabilidad_acumulada:
                seleccionados.append(matriz_puntuada[i][0])  # Agregar el individuo seleccionado
                break
    
    return seleccionados

def seleccion_por_torneo(matriz_puntuada, num_seleccionados, tam_torneo):
    seleccionados = []
    # Repetir el torneo para seleccionar el número deseado de individuos
    for _ in range(num_seleccionados):
        # Seleccionar aleatoriamente un subconjunto de individuos para el torneo
        participantes = random.sample(matriz_puntuada, tam_torneo)
        # Seleccionar al mejor individuo dentro del torneo
        mejor_individuo = max(participantes, key=lambda x: x[1])
        seleccionados.append(mejor_individuo[0])  # Agregar el individuo seleccionado
        
    return seleccionados

    
def partir_cromosomas(individuo): ## Todos los bits destinados al combustible se segmentan en un único cromosoma
    cromosoma1 = individuo[:4]
    cromosoma2 = individuo[4:8]
    cromosoma3 = individuo[8:12]
    cromosoma4 = individuo[12:15]
    cromosoma5 = individuo[15:]
    return cromosoma1, cromosoma2, cromosoma3, cromosoma4, cromosoma5

todos_los_hijos = []

def cruzamiento(seleccionados, num_descendencia, matriz_ordenada):
    global todos_los_hijos
    hijos_cruce = []
    candidatos = seleccionados.copy()
    matriz_idx = len(seleccionados)

    while num_descendencia > 0:
        if len(candidatos) < 2:
            break
        
        padre1 = random.choice(candidatos)
        padre2 = random.choice(candidatos)
        
        # Asegurarse de que padre1 y padre2 sean diferentes
        while padre1 == padre2:
            padre2 = random.choice(candidatos)
        
        # Realizar cruzamiento
        hijo = realizar_cruzamiento2(padre1, padre2, matriz_ordenada, matriz_idx)
        
        if hijo not in todos_los_hijos:
            hijos_cruce.append(hijo)
            todos_los_hijos.append(hijo)
            num_descendencia -= 1  # Decrementar el contador de descendencia

    return hijos_cruce

def realizar_cruzamiento2(padre1, padre2, matriz_ordenada, matriz_idx, max_intentos=100):
    global todos_los_hijos, historico
    cromosoma11, cromosoma21, cromosoma31, cromosoma41, cromosoma51 = partir_cromosomas(padre1)
    cromosoma12, cromosoma22, cromosoma32, cromosoma42, cromosoma52 = partir_cromosomas(padre2)

    intentos = 0
    while intentos < max_intentos:
        cromosoma13 = random.choice([cromosoma11, cromosoma12])
        cromosoma23 = random.choice([cromosoma21, cromosoma22])
        cromosoma33 = random.choice([cromosoma31, cromosoma32])
        cromosoma43 = random.choice([cromosoma41, cromosoma42])
        cromosoma53 = random.choice([cromosoma51, cromosoma52])
        hijo = cromosoma13 + cromosoma23 + cromosoma33 + cromosoma43 + cromosoma53
        hijo_tupla = tuple(hijo)

        # Verificar si el hijo es diferente a los padres y no está en la lista global de hijos ni en el historico
        if hijo != padre1 and hijo != padre2 and hijo not in todos_los_hijos and hijo_tupla not in historico:
            return hijo
        
        intentos += 1

    # Si no se puede generar un nuevo hijo después del número máximo de intentos, usar siguiente mejor individuo
    while matriz_idx < len(matriz_ordenada):
        siguiente_individuo = matriz_ordenada[matriz_idx][0]
        matriz_idx += 1
        hijo = realizar_cruzamiento2(padre1, siguiente_individuo, matriz_ordenada, matriz_idx, max_intentos)
        if hijo not in todos_los_hijos:
            return hijo

    # Si no se puede generar un nuevo hijo devuelve uno de los padres
    return random.choice([padre1, padre2])


def comprobar_individuos(individuo):
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
    global todos_los_individuos_mutados
    individuos_mutados = []
    individuos_a_mutar = random.sample(individuos, n_a_mutar)  # Seleccionar al azar n_a_mutar individuos
                                                                
    for individuo in individuos_a_mutar:
        intentos = 0
        max_intentos = 200  # Limitar el número de intentos de mutación
        while intentos < max_intentos:
            individuo_mutado = list(individuo)  # Convertir el individuo a una lista mutable
            for i in range(len(individuo_mutado) - 18):
                if random.random() < tasa_mutacion:  # Comprobar si se debe mutar este bit
                    individuo_mutado[i] = 1 - individuo_mutado[i]  # Mutación: cambiar 0 a 1 o 1 a 0
            
            if comprobar_individuos(individuo_mutado) and individuo_mutado not in todos_los_individuos_mutados:
                individuos_mutados.append(individuo_mutado)
                todos_los_individuos_mutados.append(individuo_mutado)
                break
            intentos += 1
        
        if intentos == max_intentos:
            nuevo_individuo = random.choice(individuos)
            # Realizar mutación en el nuevo individuo seleccionado
            nuevo_individuo_mutado = mutacion([nuevo_individuo], tasa_mutacion, 1)[0]
            individuos_mutados.append(nuevo_individuo_mutado)
            todos_los_individuos_mutados.append(nuevo_individuo_mutado)
    
    return individuos_mutados

def modificar_poblacion(poblacion, seleccionados, hijos_cruce, individuos_mutados, nuevo_individuo):
    # Mantener los dos individuos seleccionados originalmente
    nueva_poblacion = list(seleccionados)
    
    # Agregar los hijos resultantes del cruce
    nueva_poblacion.extend(hijos_cruce)
    
    # Agregar los individuos mutados
    nueva_poblacion.extend(individuos_mutados)
    
    # Agregar los individuos creados para completar la población
    nueva_poblacion.extend(nuevo_individuo)
    
    for individuo in seleccionados:
        historial_origen.append('Selección')
        
    # Actualizar el historial de origen para los hijos resultantes del cruce
    for _ in hijos_cruce:
        historial_origen.append('Cruce')
   
    # Actualizar el historial de origen para los individuos mutados
    for _ in individuos_mutados:
        historial_origen.append('Mutación')
    
    # Actualizar el historial de origen para el nuevo individuo
    for _ in nuevo_individuo:
        historial_origen.append('Creación')
       
    
    return nueva_poblacion


maximo_generaciones=10

# Inicializar el contador de generaciones
numero_de_generaciones = 0
num_individuos =10


num_seleccionados=4
num_descendencia=2
n_a_mutar0=2
tasa_mutacion=0.15
num_individuos_a_crear=2

poblacion = []
historico={}
individuos = creacion_individuos(num_individuos, 33, 6)
historial_origen = ['Creación' for _ in range(num_individuos)]
poblacion += individuos
historico_poblacion=[]   
while numero_de_generaciones < maximo_generaciones:
    n_a_mutar=n_a_mutar0
    print("Generación:", numero_de_generaciones)
    matriz_puntuada = Evaluacion_individuos(poblacion)
    matriz_ordenada = ordenar_matriz_puntuacion(matriz_puntuada)
    for tupla in matriz_ordenada:
        individuo = tupla[0]
        puntuacion = tupla[1]
        Q=tupla[2]
        W=tupla[3]
        CO_emission=tupla[4]
        CO2_emission=tupla[5]
        individuo_tupla = tuple(individuo)  # Convertimos la lista en una tupla, ya que da error al acceder a los términos
        if individuo_tupla not in historico: 
            historico[individuo_tupla] = (puntuacion, numero_de_generaciones ,Q, W, CO_emission, CO2_emission)
    historico_poblacion.append((matriz_puntuada))   
    seleccionados = seleccion(matriz_ordenada, num_seleccionados) 
    hijos_cruce = cruzamiento(seleccionados, num_descendencia, matriz_ordenada)
    if len(hijos_cruce) < num_descendencia:
        n_a_mutar+= 1
    mutados = mutacion(seleccionados, tasa_mutacion, n_a_mutar)
    nuevo_individuo = creacion_individuos(num_individuos_a_crear, 33, 6)
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

    # Incrementar el contador de generaciones
    numero_de_generaciones += 1
    
    # Para evitar fallos al guardar los archivos pongo la fecha con hora en el nombre 
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")  # YYYYMMDD_HHMMSS

    # Verificar la condición de parada 
    if numero_de_generaciones >= maximo_generaciones:
        columnas = ["Individuos", "Fitness", "Generación", "Calor Liberado", "Trabajo de Expansión", "Emisiones de CO", "Emisiones de CO2"]
    
        
        now = datetime.datetime.now()
        timestamp = now.strftime("%Y%m%d_%H%M%S")  # YYYYMMDD_HHMMSS
        
        # Nombres de los archivos CSV
        nombre_prueba_csv = f"Historico_{timestamp}.csv"
        nombre_origen_csv = f"Origen_{timestamp}.csv"
        
        # Exportar datos de la variable historico
        with open(nombre_prueba_csv, "w", newline='') as archivo_csv:
            escritor_csv = csv.writer(archivo_csv)
            archivo_csv.write(f"# Máximo de generaciones: {maximo_generaciones}, Número de generaciones: {numero_de_generaciones}, Número de individuos: {num_individuos},\n")
            archivo_csv.write(f"# Número de seleccionados: {num_seleccionados}, Número de descendencia: {num_descendencia}, Número a mutar: {n_a_mutar}, Tasa de mutación: {tasa_mutacion}, Número de individuos a crear: {num_individuos_a_crear}\n")
    
            
            # Escribir los encabezados de las columnas
            escritor_csv.writerow(columnas)
            
            # Escribir los datos de la matriz ordenada
            for generacion_idx, generacion in enumerate(historico_poblacion):
                for individuo in generacion:
                    cromosoma, fitness, Q, W, CO_emission, CO2_emission = individuo
                    fila = [cromosoma, fitness, generacion_idx, Q, W, CO_emission, CO2_emission]
                    escritor_csv.writerow(fila)
                
                # Exportar datos  de historial_origen (Origen de los individuos)
        with open(nombre_origen_csv, "w", newline='') as archivo_csv:
            escritor_csv = csv.writer(archivo_csv)
            
           
            escritor_csv.writerow(["Historial Origen"])
            
            # Escribe los datos de historial_origen 
            for valor in historial_origen:
                escritor_csv.writerow([valor])
