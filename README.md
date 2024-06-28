# Sudoku Genetic Algorithm

## Integrantes
- [Daysi Lizet Campos Fernández](https://github.com/BangBanguwu)
- [Jesus Mauricio Huayhua Flores](https://github.com/JesusHuayhua)
- [Juan Adolfo Mendiz Gamarra](https://github.com/Medz-10)
- [Flavio Roberto Pujay Angeles](https://github.com/JAYPU-11)

## Resumen

Este trabajo aborda la resolución de un rompecabezas Sudoku utilizando un algoritmo genético.
El objetivo general es encontrar una solución válida para el Sudoku.
Se planteó una solución mediante la implementación de un algoritmo genético que incluye clases
para representar candidatos, manejar poblaciones, realizar cruces y mutaciones, y seleccionar
los mejores individuos mediante torneos. El algoritmo evoluciona la población de soluciones potenciales
a través de la selección, el cruce y la mutación de candidatos hasta encontrar una solución válida
o alcanzar el número máximo de generaciones. Los resultados muestran que el algoritmo puede encontrar
soluciones válidas para el Sudoku, demostrando su efectividad. La conclusión principal es que el algoritmo
genético es una herramienta viable y eficiente para resolver problemas de Sudoku.


## Compilacion

```bash
g++ main.cpp -o main.exe
./main.exe > salida.txt
```