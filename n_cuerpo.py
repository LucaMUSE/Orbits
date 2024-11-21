from numpy import array , matrix, reshape, zeros

def n_cuerpos (U,t):
    
    #Nc = dimension en el espacio
    #Nb = numero de body
    #2 posicion = 0; velocidad = 1
    
    pv = reshape (U,(Nb,Nc,2)) 
    r = reshape (U[:,:,0], (Nb,Nc)) 
    v = reshape (U[:,:,1], (Nb,Nc)) 
    
    Fs = zeros(2*Nb*Nc)
    pFs = reshape (Fs,(Nb,Nc,2))
    drdt = reshape(pFs[:,:,0],(Nb,Nc)) 
    dvdt = reshape(pFs[:,:,1],(Nb,Nc))
    drdt = v 
    

