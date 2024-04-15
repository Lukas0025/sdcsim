#!/usr/bin/python3

from copy import deepcopy
import subprocess
from random import randint

#SIMD DNA PARAMETRY
DOMAINS       = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
CELL_SIZE     = 3
INTS_COUNT    = 30
NUM_OF_STATES = 2

INITALS       = [
    "000",
    "001",
    "010",
    "011",
    "100",
    "101",
    "110",
    "111"
]

FINALS        = [
    "001",
    "010",
    "011",
    "100",
    "101",
    "110",
    "111",
    "111",
]

INITAL_STATES = None

# evoluce
POPSIZE   = 8
TOUR = 2            # baze turnajove selekce GA - doporucuji nemenit
PMUT = 4            # pravdepodobnost mutace kazdeho jednotliveho genu <rule>
# PCROSS = 70       # nyni se krizeni nevyuziva (neimplementovano)
MGENES = 2          # pocet mutovanych hdnot na chromozom
GENERATIONS = 1000  # max. pocet generaci GA



def generateCell(strands):
    
    cellString = ""
    isLower    = False
    pos        = 0

    while (pos < CELL_SIZE):
        isUpper    = False
        for strand in strands:
            if len(strand) >= 2 and pos < CELL_SIZE - 2:
                if strand[0] == DOMAINS[pos] + "*" and strand[0] == DOMAINS[pos] + "*":
                    #bind it
                    if isLower:
                        cellString += "}"

                    cellString += "["
                    for i in range(len(strand)):
                        if strand[i] == DOMAINS[pos] + "*":
                            cellString += DOMAINS[pos]
                        else:
                            break

                        pos += 1
                        if pos >= CELL_SIZE:
                            break

                    cellString += "]"

                    isUpper    = True
                    isLower    = False

                    break

        if not(isUpper):
            if not(isLower):
                cellString += "{"

            cellString += DOMAINS[pos]

            pos     += 1
            isLower = True

    
    if isLower:
        cellString += "}"

    return cellString
                             

def mutate(strand):
    method = randint(0, 2) # 0 - insertion, 1 - deletition, 3 - substitution

    if len(strand) == 0:
        domain     = randint(0, len(DOMAINS) - 1)
        complement = randint(0, 1)

        if complement:
            complement = "*"
        else:
            complement = ""

        strand.append(DOMAINS[domain] + complement)

        return strand

    if method == 0:
        pos        = randint(0, len(strand))
        domain     = randint(0, len(DOMAINS) - 1)
        complement = randint(0, 1)

        if complement:
            complement = "*"
        else:
            complement = ""

        strand.insert(pos, DOMAINS[domain] + complement)

    elif method == 1:
        pos = randint(0, len(strand) - 1)
        strand.pop(pos)

    elif method == 2:
        pos        = randint(0, len(strand) - 1)
        domain     = randint(0, len(DOMAINS) - 1)
        complement = randint(0, 1)

        if complement:
            complement = "*"
        else:
            complement = ""

        strand[pos] = DOMAINS[domain] + complement

    return strand

def complement(ch):
    if ch[-1] == "*":
        return ch[0]
    
    return ch + "*"

class Genome():
    def __init__(self, number_of_states, number_of_instructions, states = None):
        if states is None:
            self.datas   = [[[]] for _ in range(number_of_states)]
        else:
            self.datas   = states[:]
        self.insts   = [[[]] for _ in range(number_of_instructions)]
        self.fitness = 0

    def getMacros(self):
        data = ""
        for i in range(len(self.datas)):
            data += "    " + str(i) + " " + generateCell(self.datas[i]) + "\n"

        return data
    
    def getInsts(self):
        data = ""
        for i in range(len(self.insts)):
            data += "    " + "".join("{" + "".join(v) + "} " for v in self.insts[i]) + " # instruction " + str(i) + "\n"

        return data

    def getAsm(self):
        code =  "define:\n"
        code += self.getMacros()
        code += "\ndata:\n    "
        code += "\n    ".join(INITALS)
        code += "\n\ninstructions:\n"
        code += self.getInsts()

        return code

    def toUpperDomains(self, strand):
        upper = []

        isLower = False
        isUpper = False
        isDs    = False

        for ch in strand:
            if ch == "[":
                isDs = True
            elif ch == "]":
                isDs = False
            elif ch == "{":
                isLower = True
            elif ch == "}":
                isLower = False
            elif ch == ".":
                continue
            elif ch == "<":
                isUpper = True
            elif ch == ">":
                isUpper = False
            elif not(isUpper):
                if ch.isnumeric(): # is macro
                    upper += self.toUpperDomains(generateCell(self.datas[int(ch)]))
                elif ch == "*" and not(isLower):
                    upper[-1] += "*"
                else:
                    if isLower:
                        upper.append("-")
                    else:
                        upper.append(complement(ch))

        return upper

    def calcFitness(self, out, ref):
        self.fitness = 0

        for i in range(len(out)):
            tOut = self.toUpperDomains(out[i])
            tRef = self.toUpperDomains(ref[i])

            for j in range(len(tOut)):
                if tOut[i] == tRef[i]:
                    self.fitness += 1

        for i in range(len(self.datas)):
            for j in range(len(self.datas)):
                if i != j:
                    if self.toUpperDomains(generateCell(self.datas[i])) != self.toUpperDomains(generateCell(self.datas[j])):
                        self.fitness *= 2
                    else:
                        self.fitness /= 2

    def mutate(self):
        
        #
        # First mutate datas
        #

        method = randint(0, 2)
        state  = randint(0, len(self.datas) - 1)

        if method == 0:
            self.datas[state].append(mutate([]))
        elif method == 1:
            if len(self.datas[state]) != 0:
                index = randint(0, len(self.datas[state]) - 1)
                self.datas[state].pop(index)

        elif method == 2:
            if len(self.datas[state]) != 0:
                index = randint(0, len(self.datas[state]) - 1)
                self.datas[state][index] = mutate(self.datas[state][index])

        #
        # Second mutate instructions
        #

        method = randint(0, 2)
        ints   = randint(0, len(self.insts) - 1)

        if method == 0:
            self.insts[ints].append(mutate([]))
        elif method == 1:
            if len(self.insts[ints]) != 0:
                index = randint(0, len(self.insts[ints]) - 1)
                self.insts[ints].pop(index)
        elif method == 2:
            if len(self.insts[ints]) != 0:
                index = randint(0, len(self.insts[ints]) - 1)
                self.insts[ints][index] = mutate(self.insts[ints][index])

def eval_pop(pop):
    best = pop[0]

    for ind in pop:
        asm = ind.getAsm()
        
        f = open("ind.sdcams", "w")
        f.write(asm)
        f.close()

        out = subprocess.check_output(['./build/sdcsim', 'ind.sdcams', '-d', '--silent'])
        out = out.decode("ascii").split("\n")
        
        tOuts = []

        for i in range(len(INITALS)):
            tOuts.append(out[i - len(INITALS) - 1])

        out = tOuts

        ind.calcFitness(out, FINALS)

        if (ind.fitness > best.fitness):
            bets = ind
    
    return best



def GA(max_gen):
    
    generation = 0
    stop_flag  = False
    total_best = 0

    pool1 = [Genome(NUM_OF_STATES, INTS_COUNT, INITAL_STATES) for _ in range(POPSIZE)]
    pool2 = [Genome(NUM_OF_STATES, INTS_COUNT, INITAL_STATES) for _ in range(POPSIZE)]

    while not(stop_flag) and generation < max_gen:

        generation = generation + 1
        
        if generation & 1:
            population = pool1
            next_population = pool2
        else:
            population = pool2
            next_population = pool1

        best = eval_pop(population)

        if total_best < best.fitness:
            print(f"GEN {generation} Best fitness: {best.fitness}")
            total_best = best.fitness

        # elitizmus
        next_population[0] = deepcopy(best)
        next_population[1] = deepcopy(best)

        for _ in range(MGENES):
            next_population[1].mutate()

        # tvorba nove populace (jedinci na pozicich [0] a [1] jsou jiz
        # vybrani v ramci elistismu vyse)
        for i in range(2, POPSIZE):
            next_population[i].fitness = 0 # vitez jde primo do nove populace
            for t in range(0, TOUR):
                c = randint(0, POPSIZE - 1)
                if population[c].fitness >= next_population[i].fitness:
                    next_population[i] = deepcopy(population[c])

            # mutace
            for _ in range(MGENES):
                next_population[i].mutate()

    return stop_flag


GA(GENERATIONS)