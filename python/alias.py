from numpy import array , matrix, reshape

#alias en vectores

V = array ([1,2,3])
pV = V #alias
pV [0] = 10

print (V)

#sale el mismo numero por eso son alias
print (id(V))
print (id(pV))

U = V.copy() #CLONING
U [0] = 18 #si hago el copy y voy a cambiar el valor del copy no noy a cambiar tambien el valor de el elemento copiado
print (V)

U = array ([1,2,3,4])
pU = reshape (U,(2,2))
pU [0,0] = 8
print(pU)

#habemos defindo punteros y punteros de punteros
