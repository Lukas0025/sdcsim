#!/usr/bin/python3

from copy import deepcopy
import subprocess
from random import randint
import argparse
import os
import tempfile

parser = argparse.ArgumentParser(
                    prog='evo.py',
                    description='Simple evolutionary SIMD||DNA designer',
                    epilog='Created by Lukas Plevac <xpleva07> BUT FIT')

parser.add_argument('--pop',                     type=int, help='Size of population') 
parser.add_argument('--pmut',                    type=int, help='Probability of mutation for one gene') 
parser.add_argument('--tour',                    type=int, help='Baze turnajove selekce GA') 
parser.add_argument('--mgenes',                  type=int, help='Number of mutated gene values') 
parser.add_argument('--gen',                     type=int, help='Number of generations') 
parser.add_argument('--drand',                   type=int, help='Randomize data every N generations if set') 
parser.add_argument('--profile',      action='store_true', help='Do profile with output csv') 

args = parser.parse_args()

PROFILE = args.profile

#SIMD DNA PARAMETRY
DOMAINS       = ["A", "B", "C", "D", "E", "F", "G", None]

CELL_LEN      = 4

STRANDS_COUNT = 2
STRANDS_LEN   = 4
INTS_COUNT    = 3

NUM_OF_STATES = 2

INITALS       = [
    "00",
    "01",
    "10",
    "11"
]

FINALS        = [
    "00",
    "00",
    "10",
    "10"
]

STATES = None

# evoluce
POPSIZE        = args.pop     # Velikost populace
TOUR           = args.tour    # baze turnajove selekce GA
PMUT           = args.pmut    # pravdepodobnost mutace kazdeho jednotliveho genu strand or data
MGENES         = args.mgenes  # pocet mutovanych hodnot v genu
GENERATIONS    = args.gen     # max. pocet generaci GA
RANDOMIZE_DATA = args.drand   # Randomize states reprezentation every N generations

#compute max fintess
MAX_FITNESS = len(FINALS) * ((len(FINALS[0]) * CELL_LEN) ** 2)

def getInsStrand(gen):
    offset     = gen[0]
    length     = gen[1]
    end        = gen[2]
    complement = gen[3]
    use        = gen[4]

    strand = []

    if not(use):
        return strand

    for i in range(length):
        strand.append(
            DOMAINS[(i + offset) % CELL_LEN]
        )

        if complement:
            strand[-1] += "*"

    if end is not None:
        strand.append(end)

        if complement:
            strand[-1] += "*"

    return strand


def complement(ch):
    if ch[-1] == "*":
        return ch[0]
    
    return ch + "*"

def randomDomain():
    isComplement = randint(0, 1)
    domain       = randint(0, len(DOMAINS) - 1)

    if isComplement and not(DOMAINS[domain] is None):
        return complement(DOMAINS[domain])

    return DOMAINS[domain]

def randomBinding():
    return randint(0, 2) # 0 - not binded, 1 - binded, 2 - binded but end of current strand

def randomStrand(MaxLen):
    strand = [
        randint(0, CELL_LEN), # offset
        randint(0, MaxLen), # len
        DOMAINS[randint(0, len(DOMAINS) - 1)], # end
        randint(0, 1), # complement
        randint(0, 1) # use
    ]

    return strand

def randomStrands(count, length):
    strands = []
    for _ in range(count):
        strands.append(randomStrand(length))

    return strands

def randomData(length):
    strand = []
    for _ in range(length):
        strand.append(randomBinding())

    return strand

def toUpperBind(strand):
    upper = []

    isLower = False
    isUpper = False
    isDs    = False

    for ch in strand:
        if ch == "[":
            isDs = True
        elif ch == "]":
            isDs = False
            if len(upper) > 0:
                upper[-1] = 2
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
        elif ch == "*":
            continue
        elif not(isUpper):
            if isLower:
                upper.append(0)
            elif isDs:
                upper.append(1)

    return upper


class Genome():
    def __init__(self, number_of_states, len_of_state, number_of_instructions, number_of_insts_strands, len_of_inst_strand, states = None):
        self.datas   = [randomData(len_of_state) for _ in range(number_of_states)]
        self.mutData = True

        if states is not None:
            self.datas = states[:]
            self.mutData = False


        self.number_of_states = number_of_states
        self.len_of_state     = len_of_state

        self.insts   = [randomStrands(number_of_insts_strands, len_of_inst_strand) for _ in range(number_of_instructions)]
        self.fitness = 0

    def randomizeData(self):
        self.datas = [randomData(self.len_of_state) for _ in range(self.number_of_states)]

    def getData(self, i, alfa):

        out = []
        isUpper = False
        isLower = False
        key = 0

        for item in self.datas[i]:
            if item == 1:
                if isLower:
                    out.append("}")
                if not(isUpper):
                    out.append("[")
                isUpper = True
                isLower = False
                
                out.append(alfa[key])
            elif item == 2:
                if isLower:
                    out.append("}")
                if not(isUpper):
                    out.append("[")
                isUpper = False
                isLower = False
				
                out.append(alfa[key])
                out.append("]")
                
            else:
                if isUpper:
                    out.append("]")
                if not(isLower):
                    out.append("{")
                isUpper = False
                isLower = True
                
                out.append(alfa[key])

            key += 1

        if isUpper:
            out.append("]")
        
        if isLower:
            out.append("}")
        
        return "".join(out)
    
    def getStableData(self, i):
        out = self.datas[i][:]

        for j in range(len(out)):
            if (j + 1 < len(out)) and (out[j+1] == 2 or out[j+1] == 1) and out[j] != 2:
                continue

            if (j - 1 >= 0) and out[j-1] == 1:
                continue

            out[j] = 0


        return out

    def getMacros(self, alfa):
        data = ""
        for i in range(len(self.datas)):
            data += "    " + str(i) + " " + self.getData(i, alfa) + "\n"

        return data
        
    def getInst(self, i):
        out = ""
        for strand in self.insts[i]:
            out += "{"
            for domain in getInsStrand(strand):
                if domain is None:
                    break
                    
                out += domain	
            out += "} "
        
        return out
    
    def getInsts(self):
        data = ""
        for i in range(len(self.insts)):
            data += "    " + self.getInst(i) + " # instruction " + str(i) + "\n"

        return data

    def getAsm(self, alfa):
        code =  "define:\n"
        code += self.getMacros(alfa)
        code += "\ndata:\n    "
        code += "\n    ".join(INITALS)
        code += "\n\ninstructions:\n"
        code += self.getInsts()

        return code

    def asmXor(self, strandA, strandB, alfa): # returns number of same strands
        
        # first convert to asm mithout macros
        for i in range(len(self.datas)):
            strandA = strandA.replace(str(i), self.getData(i, alfa))
            strandB = strandB.replace(str(i), self.getData(i, alfa))

        upperA = toUpperBind(strandA)
        upperB = toUpperBind(strandB)
        pos    = 0
        score  = 0

        while True:

            if pos + CELL_LEN <= len(upperA):
                cellSame = True
                for i in range(CELL_LEN):
                    if upperA[pos + i] != upperB[pos + i]:
                        cellSame = False
                        break

                if cellSame:
                    score += CELL_LEN
            else:
                break

            pos += CELL_LEN

        return score
    
    def asmDiff(self, strandA, strandB, alfa): # returns number of same strands
        
        # first convert to asm mithout macros
        for i in range(len(self.datas)):
            strandA = strandA.replace(str(i), self.getData(i, alfa))
            strandB = strandB.replace(str(i), self.getData(i, alfa))

        upperA = toUpperBind(strandA)
        upperB = toUpperBind(strandB)

        # compute overhang penaly
        isUpper = False
        penalyA = 0

        for ch in strandA:
            if ch == "<":
                isUpper = True
            elif ch == ">":
                isUpper = False
            elif isUpper and ch != "*":
                penalyA += 1

        penalyB = 0

        for ch in strandB:
            if ch == "<":
                isUpper = True
            elif ch == ">":
                isUpper = False
            elif isUpper and ch != "*":
                penalyB += 1

        score = -abs(penalyA - penalyB)

        for i in range(len(upperA)):
            if upperA[i] == upperB[i]:
                score += 1

        return score

    def dataStability(self):
        stability = 0

        for i in range(len(self.datas)):
            local  = 0
            strand = self.getStableData(i)

            for j in range(len(strand)):
                if strand[j] == self.datas[i][j]:
                    local += 1

            stability += local / len(strand)

        return stability / len(self.datas)

    def calcFitness(self, out, ref, alfa):
        self.fitness = 0

        for i in range(len(out)):
            self.fitness += max(self.asmDiff(out[i], ref[i], alfa), 0) * self.asmXor(out[i], ref[i], alfa) * self.dataStability()

        for i in range(len(self.datas)):
            for j in range(len(self.datas)):
                if i != j:
                    a = ",".join(str(v) for v in self.getStableData(i)).replace("2", "1")
                    b = ",".join(str(v) for v in self.getStableData(j)).replace("2", "1")
                    if a == b:
                        self.fitness = 0
                        break


    def mutate(self, pmut = 4, mutcout = 2):
        if self.mutData:
            for i in range(len(self.datas)):
                if randint(0, 100) <= pmut:
                    for _ in range(mutcout):
                        pos = randint(0, len(self.datas[i]) - 1)
                        self.datas[i][pos] = randomBinding()

        for i in range(len(self.insts)):
            for j in range(len(self.insts[i])):
                if randint(0, 100) <= pmut:
                    for _ in range(mutcout):
                        pos = randint(0, 4)
                        if pos == 0:
                            self.insts[i][j][0] += randint(-1, 1) # offset
                        elif pos == 1:
                            self.insts[i][j][1] += randint(-1, 1) # length
                            self.insts[i][j][1] = max(self.insts[i][j][1], 0)
                            self.insts[i][j][1] = min(self.insts[i][j][1], STRANDS_LEN)
                        elif pos == 2:
                            self.insts[i][j][2] = DOMAINS[randint(0, len(DOMAINS) - 1)] # end
                        elif pos == 3:
                            self.insts[i][j][3] = not(self.insts[i][j][3]) # complement
                        else:
                            self.insts[i][j][4] = not(self.insts[i][j][4]) # use


def eval_pop(pop):
    best = pop[0]

    for ind in pop:
        asm = ind.getAsm(DOMAINS)
        
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(asm.encode("ascii"))
        f.close()

        out = subprocess.check_output(['./build/sdcsim', f.name, '-d', '--silent', '-t', '20'])
        out = out.decode("ascii").split("\n")
            
        tOuts = []

        for i in range(len(INITALS)):
            tOuts.append(out[i - len(INITALS) - 1])

        out = tOuts

        ind.calcFitness(out, FINALS, DOMAINS)

        if ind.fitness > best.fitness:
            best = ind

        os.unlink(f.name)

    return best

def GA(max_gen):

    global cur_insts
    
    generation = 0
    stop_flag  = False
    best       = None
    lastFit    = -1

    pool1 = [Genome(NUM_OF_STATES, CELL_LEN, INTS_COUNT, STRANDS_COUNT, STRANDS_LEN, STATES) for _ in range(POPSIZE)]
    pool2 = [Genome(NUM_OF_STATES, CELL_LEN, INTS_COUNT, STRANDS_COUNT, STRANDS_LEN, STATES) for _ in range(POPSIZE)]

    while not(stop_flag) and generation < max_gen:

        generation = generation + 1
        
        if generation & 1:
            population = pool1
            next_population = pool2
        else:
            population = pool2
            next_population = pool1

        best = eval_pop(population)

        if best.fitness > lastFit or generation % 100 == 0:
            print(f"GEN {generation} Best fitness: {best.fitness} / {MAX_FITNESS}")

        if best.fitness > lastFit:
            lastFit = best.fitness

        if best.fitness == MAX_FITNESS:
            stop_flag = True
            break

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
            next_population[i].mutate(PMUT, MGENES)

            if RANDOMIZE_DATA is not None and generation <= RANDOMIZE_DATA:
                next_population[i].randomizeData()

    return best, stop_flag

history = []
for _ in range(50):
    best, sf = GA(GENERATIONS)

    history.append([
        POPSIZE,
        TOUR,
        PMUT,
        MGENES,
        GENERATIONS,
        RANDOMIZE_DATA,
        best.fitness
    ])

    if sf:
        print("Sdcasm is:")
        print()
        print(best.getAsm(DOMAINS))

        f = open("solution.sdcams", "w")
        f.write(best.getAsm(DOMAINS))
        f.close()

    if sf and not(PROFILE):
        break

if PROFILE:
    with open("evo_output.csv", "a+") as txt_file:
        for line in history:
            txt_file.write(";".join([str(lin) for lin in line]) + "\n")