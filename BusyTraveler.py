from __future__ import print_function
from ortools.linear_solver import pywraplp
from landmarksHelper import recover_distances, recover_landmarks
import random
from time import time

# impongo il massimo tempo di computazione a 60 secondi
MAX_COMPUTING_TIME = 60 * 60
NUM_THREADS = 4

# Policy
ALL = 0
MAX_POPOLARITA = 1
MIN_DISTANZA = 2
MIXED = 3
n = 0
V = []
popolarita = []
monumenti = []
D = []
policy = "popolarita"


# Funzione ch eimplementa la parte non ricorsiva del quick sort
def quad_partition(vect1, vect2, vect3, vect4, low, high):
    i = (low - 1)
    pivot = vect1[high]

    for j in range(low, high):

        if vect1[j] <= pivot:
            i = i + 1
            vect1[i], vect1[j] = vect1[j], vect1[i]
            vect2[i], vect2[j] = vect2[j], vect2[i]
            vect3[i], vect3[j] = vect3[j], vect3[i]
            vect4[i], vect4[j] = vect4[j], vect4[i]

    vect1[i + 1], vect1[high] = vect1[high], vect1[i + 1]
    vect2[i + 1], vect2[high] = vect2[high], vect2[i + 1]
    vect3[i + 1], vect3[high] = vect3[high], vect3[i + 1]
    vect4[i + 1], vect4[high] = vect4[high], vect4[i + 1]

    return (i + 1)


# Funzione che ordina entrambi i vettori basandosi sui valori del primo
# Quick sort
def quad_sort_rec(vect1, vect2, vect3, vect4, low, high):
    if low < high:
        # pi is partitioning index, arr[p] is now 
        # at right place 
        pi = quad_partition(vect1, vect2, vect3, vect4, low, high)

        # Separately sort elements before 
        # partition and after partition 
        quad_sort_rec(vect1, vect2, vect3, vect4, low, pi - 1)
        quad_sort_rec(vect1, vect2, vect3, vect4, pi + 1, high)


def quad_sort(vect1, vect2, vect3, vect4):
    quad_sort_rec(vect1, vect2, vect3, vect4, 0, len(vect1) - 1)


# Funzione ch eimplementa la parte non ricorsiva del quick sort
def partition(vect1, vect2, vect3, low, high):
    i = (low - 1)
    pivot = vect1[high]

    for j in range(low, high):

        if vect1[j] <= pivot:
            i = i + 1
            vect1[i], vect1[j] = vect1[j], vect1[i]
            vect2[i], vect2[j] = vect2[j], vect2[i]
            vect3[i], vect3[j] = vect3[j], vect3[i]

    vect1[i + 1], vect1[high] = vect1[high], vect1[i + 1]
    vect2[i + 1], vect2[high] = vect2[high], vect2[i + 1]
    vect3[i + 1], vect3[high] = vect3[high], vect3[i + 1]
    return (i + 1)


# Funzione che ordina entrambi i vettori basandosi sui valori del primo
# Quick sort
def triple_sort_rec(vect1, vect2, vect3, low, high):
    if low < high:
        # pi is partitioning index, arr[p] is now 
        # at right place 
        pi = partition(vect1, vect2, vect3, low, high)

        # Separately sort elements before 
        # partition and after partition 
        triple_sort_rec(vect1, vect2, vect3, low, pi - 1)
        triple_sort_rec(vect1, vect2, vect3, pi + 1, high)


def triple_sort(vect1, vect2, vect3):
    triple_sort_rec(vect1, vect2, vect3, 0, len(vect1) - 1)


def set_up(pol, numMon=50):
    global n, V, popolarita, monumenti, D, policy

    # nome dei monumenti da visitare
    monumenti = recover_landmarks()

    random.seed(a="paperino", version=2)

    # vettore di popolarità dei monumenti del grafo (il primo elemento
    # è la popolarità del nodo utente, pari a zero. Inserito per semplificare
    # la gestione degli indici nell'implementazione)
    popolarita = [random.randint(0, 100) for i in range(0, len(monumenti))]

    # vettore di distanze dall'utente ai nodi
    # (da incorporare nella matrice delle distanze)
    du = [831, 623, 700, 712, 720, 825, 593, 692, 679, 717, 592, 572, 858, 515, 174, 625, 998, 681, 792, 794, 784, 1311, 625, 978, 704, 290, 713, 712, 637, 795, 684, 1008, 571, 799, 623, 541, 770, 710, 623, 1069, 696, 787, 745, 374, 800, 788, 711, 710, 583, 293, 775, 1015, 681, 503, 672, 784, 463, 666, 713, 855, 824, 1041, 314, 785, 928, 798, 733, 1025, 590, 806, 712, 359, 668, 383, 985, 1310, 1005, 1051, 720, 571, 832, 790, 824, 728, 696, 812, 1102, 590, 727, 699, 592, 737, 669, 359, 763, 769, 828, 503, 655, 1160]

    # Prima di prendere le distanze applico una policy per semplificare il modello
    if (pol == MAX_POPOLARITA):
        print("Filtrando monumenti per massima popolarità")
        triple_sort(popolarita, monumenti, du)
        popolarita = popolarita[-numMon:]
        monumenti = monumenti[-numMon:]
        du = du[-numMon:]
    elif (pol == MIN_DISTANZA):
        print("Filtrando monumenti per minima distanza")
        triple_sort(du, monumenti, popolarita)
        popolarita = popolarita[:numMon]
        monumenti = monumenti[:numMon]
        du = du[:numMon]
    elif (pol == MIXED):
        mixed = [0] * len(monumenti)
        for i in range(0, len(du)):
            mixed[i] = du[i] - 8 * popolarita[i]
        quad_sort(mixed, popolarita, du, monumenti)
        popolarita = popolarita[:numMon]
        monumenti = monumenti[:numMon]
        du = du[:numMon]
    else:
        print("Nessuno filtro applicato")
        popolarita = popolarita[:numMon]
        monumenti = monumenti[:numMon]
        du = du[:numMon]

    # matrice delle distanze, per semplificare i test assumo in tale istanza che
    # la distanza associata all'arco Xij è uguale alla distanza associata all'arco Xji
    # (le distanza sono espresse in termini di tempo)
    distanze_info = recover_distances()

    distanze = [[0 for j in range(i, numMon - 1)] for i in range(0, numMon - 1)]

    for i in range(0, len(monumenti)):
        for j in range(i + 1, len(monumenti)):
            c = 0
            while not ((distanze_info[c]["Start"] == monumenti[i] and distanze_info[c]["End"] == monumenti[j]) or (
                    distanze_info[c]["End"] == monumenti[i] and distanze_info[c]["Start"] == monumenti[j])):
                c += 1
            distanze[i][j - i - 1] = int(distanze_info[c]["Seconds"])

    distanze.append([])

    # aggiungo distanze dall'utente ai nodi
    V = set(range(len(distanze)))
    D = [[0 if i == j
          else distanze[i][j - i - 1] if j > i
    else distanze[j][i - j - 1]
          for j in V] for i in V]

    cp = du[:]
    cp.insert(0, 0)
    D.insert(0, cp)
    print(len(cp), len(D), len(V))
    for i in range(1, len(V) + 1):
        D[i].insert(0, cp[i])

    # definisco l'insieme dei nodi e la sua cardinalità
    V = set(range(1, len(monumenti)))

    # aggiungo all'insieme dei nodi, il nodo utente
    V.add(0)

    # cardinalità dell'insieme dei nodi
    n = len(V)


# tempo massimo a disposizione dell'utente
def BusyTraveler(T):
    # definisco il modello ed il suo nome
    solver = pywraplp.Solver('BusyTraveler', pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

    # definizione delle variabili X(i-j)
    infinity = solver.infinity()
    x = []
    for i in V:
        x.append([])
        for j in V:
            x[i].append(solver.IntVar(0, 1, "X({}-{})".format(i, j)))

    # definizione variabili Y(i)
    y = [solver.IntVar(0, infinity, "Y({})".format(i)) for i in V]

    # definizione variabili H(i)
    h = [solver.IntVar(0, 1, "H({})".format(i)) for i in V]

    # Vincolo di tipo (1)
    for i in V:
        solver.Add(solver.Sum(x[i][j] * D[i][j] for j in V for i in V) <= T)

    # Vincolo di tipo (2)
    for j in V:
        solver.Add(solver.Sum(x[i][j] for i in V - {j}) <= 1)

    # Vincolo di tipo (3)
    for i in V:
        solver.Add(solver.Sum(x[i][j] for j in V - {i}) <= 1)

    # Vincolo di tipo (4)
    solver.Add(y[0] == 1)

    # Vincolo di tipo (5)
    for i in V - {0}:
        solver.Add(y[i] - solver.Sum((n + 1) * x[j][i] for j in V) <= 0)

    # Vincolo di tipo (6)
    for i in V:
        for j in V:
            solver.Add(x[i][j] - y[i] <= 0)

    # Vincolo di tipo (7)
    for i in V:
        for j in V:
            solver.Add(y[j] >= y[i] + 1 - (n + 1) * (1 - x[i][j]))

    # Vincolo di tipo (8)
    for i in V:
        solver.Add(h[i] <= y[i])

    # funzione obiettivo
    # objective = solver.Objective()
    solver.Maximize(solver.Sum(h[i] * popolarita[i] for i in V))

    # TEMPO LIMITE
    solver.set_time_limit(MAX_COMPUTING_TIME * 1000)
    # numero di thread per il solver
    solver.SetNumThreads(NUM_THREADS)

    start = time()

    # optimizing
    status = solver.Solve()

    end = time()

    # stampo la soluzione
    if status == pywraplp.Solver.OPTIMAL:
        solved = "Si"
        ub = solver.Objective().Value()
        lb = ub
        f.write('Valore della funzione obbiettivo = {}\n'.format(lb))

        mons = []
        for i in V - {0}:
            if (y[i].solution_value() != 0):
                mons.append((i, y[i].solution_value()))
        mons.sort(key=lambda x: x[1])

        f.write("Percorso trovato con t = {}:\n".format(T))
        f.write("Utente")

        for mon in mons:
            f.write("----->{}".format(monumenti[mon[0]]))
        f.write(".\n")

        f.write('Problema risolto in {} minuti\n'.format((solver.wall_time() / 1000) / 60))

    else:
        lb = solver.Objective().Value()
        ub = solver.Objective().BestBound()
        solved = "No"
        f.write('Il problema non ha soluzione ottima entro {} secondi\n'.format(MAX_COMPUTING_TIME))
        f.write('Lower Bound = {}\n'.format(lb))
        f.write('Upper Bound = {}\n'.format(ub))

    f.write("--------------------\n\n")

    # scrivo risultati sul csv
    csv_results(T, lb, ub, solved, len(monumenti), end - start)

    return solved


csv_filename = "BusyTravelerResults.csv"


def csv_results(t, lb, up, solved, numMon, time_spent):
    import os

    # if the file does not exist
    if (not os.path.exists(csv_filename)):
        f = open(csv_filename, "w")
        f.write(
            "Tempo Del Viaggiatore, Lower Bound, Upper Bound, Risolto In 1 h, Monumenti Considerati, Politica Selezione Monumenti, Tempo\n")
    else:
        f = open(csv_filename, "a")
    f.write(str(t) + ", " + str(lb) + ", " + str(up) + ", " + str(solved) + ", " + str(numMon) + ", " + str(
        policy) + ", " + str(time_spent) + "\n")
    f.close()


def self_test(p, numNodes, T):
    # @p: policy da usare
    # @T: tempo considerato
    # @numNodes: monumenti iniziali inizialmente considerati
    delta = int(numNodes / 2)
    perc = open("Percentage.txt", "a")
    perc.write("Self test staring with {} num nodes\n".format(numNodes))
    perc.close()
    while True:
        if delta < 1:
            break

        set_up(p, numNodes)
        solved = BusyTraveler(T)
        if solved == "Si":
            # lascio inalterato il delta
            numNodes += delta
            perc = open("Percentage.txt", "a")
            perc.write("Solved... re-trying with {} num nodes\n".format(numNodes))
            perc.close()
            continue
        else:
            numNodes -= delta
            delta = int(delta / 2)
            perc = open("Percentage.txt", "a")
            perc.write("NOT solved... re-trying with {} num nodes\n".format(numNodes))
            perc.close()
            continue


f = open("ResultsBusyTraveller.txt", "a")
if __name__ == "__main__":

    # inizio con tutti i nodi, poi restringo se necessario
    # dentro il metodo BusyTraveller
    hours = [5]
    policies = [MIXED]


    set_up(MIXED, 10)
    solved = BusyTraveler(10000)

    # test finale, con 12ore date al solver, vedere risultati
    #self_test(MIXED, 100, 5)



    """
    for hour in hours:
        T = 60*60*hour
        perc = open("Percentage.txt", "a")
        perc.write("Hours: {}\n".format(hour))
        perc.close()
        for p in policies:
            perc = open("Percentage.txt", "a")
            perc.write("Policy: {}\n".format(p))
            perc.close()
            for i in range(15, 40):
                perc = open("Percentage.txt", "a")
                perc.write("Nodes: {}\n".format(i))
                perc.close()
                print("Currently analizing ", i)
                s = time()
                set_up(p, i)
                solved = BusyTraveler(T)
                e = time()

                if e -s > 3600:
                    break

                if solved == "Si":
                    print("Solved for ", i)
                else:
                    print("Not solved for ", i)
    """

    perc = open("Percentage.txt", "a")
    perc.write("Completed!\n")
    perc.close()
    f.close()
