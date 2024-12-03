from itertools import permutations

# Función para encontrar todas las particiones de un número N en grupos
def particiones(n):
    resultado = []
    particion = [0] * n  # Creamos una lista para almacenar la partición
    def auxiliar(m, i):
        if m == 0:
            resultado.append(particion[:i])  # Agregamos la partición encontrada
            return
        for j in range(1, m + 1):  # Exploramos las posibilidades de partición
            particion[i] = j
            auxiliar(m - j, i + 1)
    
    auxiliar(n, 0)
    return resultado

# Función para crear una descripción de la partición
def descripcion_particion(particion, N):
    descripcion = []
    for i in range(N):
        if i < len(particion) and particion[i] > 0:
            descripcion.append(f"{particion[i]} satélites en la órbita {i + 1}")
        elif i < len(particion) and particion[i] == 0:
            descripcion.append(f"ningún satélite en la órbita {i + 1}")
        else:
            descripcion.append(f"ningún satélite en la órbita {i + 1}")
    return ", ".join(descripcion)

# Función para calcular y mostrar el número de formas y las formas mismas
def calcular_formas(N):
    # Obtenemos todas las particiones de N en grupos
    particiones_n = particiones(N)
    
    formas = set()  # Usamos un conjunto para evitar duplicados
    for particion in particiones_n:
        # Agregamos ceros para completar el número de órbitas
        particion = particion + [0] * (N - len(particion))  # Completamos las órbitas con ceros
        formas.add(tuple(particion))  # Agregamos la partición como una tupla para evitar duplicados
    
    # Imprimimos el número de formas
    print(f"Número de formas: {len(formas)}")
    
    # Imprimimos las formas mismas y su descripción
    for forma in formas:
        print(forma)
        descripcion = descripcion_particion(forma, N)
        print(f"Descripción: {descripcion}")
    
    # Creamos la matriz
    matriz = []
    for forma in formas:
        matriz.append(list(forma))
    
    return matriz

# Pedimos al usuario el número de satélites
N = int(input("¿Cuántos satélites deseas lanzar a órbita? "))

# Calculamos y mostramos las formas
matriz = calcular_formas(N)

# Para visualizar la matriz
print("\nMatriz de las formas:")
for fila in matriz:
    print(fila)
