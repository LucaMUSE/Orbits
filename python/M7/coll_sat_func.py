import itertools

def calcola_distribuzione_satelliti(satelliti):
    distribuzioni = []

    # Generiamo tutte le partizioni di 'satelliti'
    for r in range(1, satelliti + 1):  # Da 1 a N orbite
        for combinazione in itertools.combinations_with_replacement(range(1, satelliti +1 ), r):
            if sum(combinazione) == satelliti:
                distribuzioni.append(tuple(sorted(combinazione)))  # Ordina per chiarezza

    # Restituisce il numero di modi diversi di collocare i satelliti
    return len(distribuzioni), distribuzioni

def main():
    satelliti = int(input("Quanti satelliti vuoi lanciare? "))
    modi, distribuzioni = calcola_distribuzione_satelliti(satelliti)
    print(f"Ci sono {modi} modi diversi di collocare i satelliti sulle orbite.")
    print("Le distribuzioni possibili sono:")
    for distribuzione in distribuzioni:
        print(distribuzione)

if __name__ == "__main__":
    main()
