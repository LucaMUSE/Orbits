
def Leap-Frog (F,t1,t2,U):
    
    t2 - t1 = dt 
    
    if t1==0:
        
    U1 = U + (dt) * F(U)
    
    return U1

else:
    
    U2 = U1 + 2 * (dt) * F(U)
    
    U = U1
    
    return U2

#quando ho un metodo multipasso devo utilizzare come vettore di stato uno schema unipasso
    
