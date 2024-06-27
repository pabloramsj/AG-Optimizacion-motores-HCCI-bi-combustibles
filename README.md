# Optimización de Motores HCCI Bi-Combustibles mediante Algoritmos Genéticos

Este repositorio contiene el código fuente y los archivos necesarios para la optimización de motores HCCI (Homogeneous Charge Compression Ignition) utilizando algoritmos genéticos. El objetivo principal de este proyecto es mejorar la eficiencia y el trabajo neto producido, así como reducir las emisiones de los motores HCCI mediante el ajuste de diversos parámetros operativos usando técnicas de optimización inspiradas en la evolución natural.

## Contenido del Repositorio

- **AG.py**: Código fuente del algoritmo genético desarrollado.
- **AG Análisis resultados.ipynb**: Archivo Jupyter Notebook que permite realizar el análisis de los datos recopilados durante la ejecución del algoritmo.
- **heptane-mehl.yaml**: Archivo que contiene el mecanismo de cinética química empleado durante la evaluación de los individuos dentro del algoritmo.

## Requisitos

Para utilizar estos archivos, es necesario clonar el repositorio y asegurarse de tener instaladas las siguientes dependencias:

- [Cantera](https://cantera.org/): Librería empleada en el modelado del motor HCCI.
- [Pandas](https://pandas.pydata.org/): Librería para análisis y manipulación de datos.
- [Plotly](https://plotly.com/): Librería para visualización de datos.
- [NumPy](https://numpy.org/): Librería para computación numérica en Python.
- [Matplotlib](https://matplotlib.org/): Librería para creación de gráficos en Python.
- [SciPy](https://scipy.org/): Librería para cálculos científicos y técnicos.
- [time](https://docs.python.org/3/library/time.html): Módulo para funciones relacionadas con el tiempo. 
- [random](https://docs.python.org/3/library/random.html): Módulo para generación de números aleatorios.
- [csv](https://docs.python.org/3/library/csv.html): Módulo para manejo de archivos CSV.
- [datetime](https://docs.python.org/3/library/datetime.html): Módulo para manipulación de fechas y horas.

## Instalación

1. Clonar el repositorio en una máquina local:
    ```bash
    git clone https://github.com/pabloramsj/AG-Optimizacion-motores-HCCI-bi-combustibles.git
    cd AG-Optimizacion-motores-HCCI-bi-combustibles
    ```

2. Asegurarse de tener instalado Cantera. Puede instalarse siguiendo las instrucciones en la [documentación oficial de Cantera](https://cantera.org/install/index.html)

## Uso

### Ejecución del Algoritmo Genético

Para ejecutar el algoritmo genético, simplemente ejecutar el archivo `AG.py`:
```bash
python AG.py
