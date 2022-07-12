# Computacional
Repositorio para guardar todo el código relacionado con la asignatura de Física Computacional (UGR)

Para ejecutar los códigos se recomienda el uso de la terminal Julia REPL para lo que es necesario tener instalado el lenguaje [Julia](https://julialang.org/). 

Es necesario clonado el repositorio y, una vez hecho esto navegar a la carpeta donde se haya clonado y ejecutar el siguiente código en la terminal:
`>julia`
Para abrir la terminal interactiva de Julia.
```
  julia>using Pkg
  julia>Pkg.activate(".")
  julia>Pkg.instantiate()
```
e instalaremos así todas las librerías necesarias. Asi estaremos listos para ejecutar lo necesario.
## Informe del Modelo de Ising
El código se encuentra en la carpeta Tarea2

Primero es necesario ejecutar todas las definciones con:
```
  julia>include("Tarea2/IsingScripts.jl")
```
Para crear una simulación de tamaño **128** con spines iniciales aleatorios y temperatura **2.3** ejecutamos
```
  julia>r = rand(Red, 128, 2.3)
```
esto guardará la red en la variable **r**, para ejecutar la simulación con una longitud de 100 pMC ejecutamos:
```
  julia>s = bucle_simulacion(r, 100)
```
esto guardará un arreglo de redes en **s** que será la evolución temporal del sistema.
```
  julia>plotear_configuracion_anim(s)
```
generá un ".mp4" de la evolución temporal del sistema.
```
  julia>plotear_EyMvspMC([s], [:blue])
  julia>save("testIsing.png", ans)
```
creará una gráfica con la energía y la magnetización en función del tiempo y la guardará.
## Informe del péndulo doble
EL código se encuentra en la carpeta Tarea4
Primero es necesario ejecutar todas las definciones con:
```
  julia>include("Tarea2/PenduloDoblePlt.jl")
```
podemos crear una simulación con:
```
  julia>s = crear_simulacion(θ1 = 0.0, θ2 = 2.0)
```
esto creará una simulación con ángulos iniciales θ1 = 0 radianes y θ2 = 2.0 radianes, para ejecutarla:
```
  julia>loopfxts!(s, 100.0)
```
que ejecutará 100 segundos de simulación. Podremos animarla con:
```
  julia>animar_simulacion(s, "testPD.mp4")
```
que creará un vídeo y lo guardará en formato ".mp4" en la carpeta "Tarea4/videosPD".
