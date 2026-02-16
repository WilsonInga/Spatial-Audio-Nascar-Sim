# Spatial Audio NASCAR Simulation 🏎️🎧

![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-orange)
![License](https://img.shields.io/badge/License-MIT-blue)
![Status](https://img.shields.io/badge/Status-Stable-green)

Este proyecto implementa un motor de renderizado de audio 3D (Binaural) en MATLAB. Transforma una señal de audio monoaural (un auto de carreras) en una experiencia estéreo inmersiva, simulando que el vehículo gira alrededor de la cabeza del oyente con cambios de azimut y elevación.

## 📋 Características Principales

* **Síntesis de HRTF Procedural:** No requiere bases de datos externas (como CIPIC). Genera las *Head-Related Impulse Responses* calculando ITD (Diferencia de tiempo interaural) e ILD (Diferencia de nivel interaural) matemáticamente.
* **Motor de Convolución Rápida:** Implementa el método **Overlap-Add** con FFT para procesar el audio en bloques, permitiendo una simulación continua y eficiente.
* **Simulación de Elevación:** Utiliza filtros FIR dependientes de la altura para modificar el contenido espectral (agudos/graves) y simular que el sonido viene de arriba o abajo.
* **Trayectoria Dinámica:** El sonido realiza 2 vueltas completas (720°) con oscilación vertical sinusoidal.

## 🏗️ Arquitectura del Sistema

El flujo de procesamiento del script `nascarCS.m` sigue la siguiente lógica:

1.  **Pre-procesamiento:**
    * Carga del audio `nascar-race.wav`.
    * Conversión a mono y normalización.
    * Segmentación en bloques (Windowing) con ventana *Tukey*.

2.  **Generación de HRIRs (Head-Related Impulse Responses):**
    * Se crea un banco de filtros para 36 ángulos (cada 10°).
    * Se calcula el retraso entre oídos (ITD) basado en una aproximación esférica de la cabeza.
    * Se aplica atenuación por sombra acústica (ILD).

3.  **Renderizado en el Bucle Principal:**
    * **Interpolación Circular:** Para lograr un movimiento suave entre los 10° pre-calculados, se interpolan las HRIRs linealmente en tiempo real según la posición exacta del auto.
    * **Efecto de Elevación:** Se aplica un filtrado `apply_elevation_synthetic` que refuerza frecuencias altas cuando la fuente está elevada.
    * **Convolución FFT:** Se convoluciona el bloque de audio actual con la HRIR interpolada.

4.  **Reconstrucción:**
    * Suma de solapamiento (Overlap-Add) para reconstruir la señal estéreo final sin clics ni cortes.

## 🚀 Requisitos y Ejecución

### Prerrequisitos
* MATLAB (Cualquier versión reciente).
* **Signal Processing Toolbox** (Necesario para funciones como `tukeywin`, `fir1`, `fft`).

### Pasos para ejecutar

1.  Clona este repositorio:
    ```bash
    git clone [https://github.com/tu-usuario/spatial-audio-nascar-sim.git](https://github.com/tu-usuario/spatial-audio-nascar-sim.git)
    ```
2.  Asegúrate de que el archivo `nascar-race.wav` esté en la misma carpeta que el script `nascarCS.m`.
3.  Abre MATLAB y navega hasta la carpeta del proyecto.
4.  Ejecuta el script principal:
    ```matlab
    >> nascarCS
    ```
5.  **Resultado:**
    * Se abrirá una barra de progreso indicando la posición (Azimut/Elevación).
    * Al finalizar, se generarán dos gráficas (Audio Original vs. Audio 3D).
    * El audio resultante se reproducirá automáticamente.

> **Nota:** Para una correcta apreciación del efecto 3D, es **obligatorio usar audífonos**.

## 📊 Visualización

El script genera una salida gráfica donde se puede comparar la forma de onda monoaural original contra la salida binaural procesada, mostrando las diferencias de amplitud entre el canal izquierdo y derecho que crean la ilusión espacial.

## ✒️ Autor

**Wilson Inga**
* [Perfil de GitHub](https://github.com/tu-usuario)

---
*Este proyecto fue desarrollado con fines educativos e investigación en procesamiento de señales de audio.*
