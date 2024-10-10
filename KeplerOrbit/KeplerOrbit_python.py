def F(x, y, xd, yd):
    F1 = xd
    F2 = yd
    F3 = -x / (x**2 + y**2)
    F4 = -y / (x**2 + y**2)

    return F1, F2, F3, F4

xn = 1
yn = 0
xdn = 0
ynd = 1

# Richiamiamo la funzione F per ottenere i valori
F1, F2, F3, F4 = F(xn, yn, xdn, ynd)

# Calcoliamo i nuovi valori
xn1 = xn + 0.1 * F1
yn1 = yn + 0.1 * F2
xdn1 = xdn + 0.1 * F3
ydn1 = ynd + 0.1 * F4

# Stampa dei risultati
print("xn1 =", xn1)
print("yn1 =", yn1)
print("xdn1 =", xdn1)
print("ydn1 =", ydn1)