import argparse
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

parser = argparse.ArgumentParser(
                    prog='stats.py',
                    description='Simple statistical evaulator for sdcsim',
                    epilog='Created by Lukas Plevac <xpleva07> BUT FIT')

parser.add_argument('sdcsim',  help='Path to sdcsim simulator')
parser.add_argument('sdcasm',  help='Path to simulated adcsm file')
parser.add_argument('-T', '--temperarure',       type=float, help='Center of tested temperatures') 
parser.add_argument('-t', '--time',              type=int, help='Center of tested times') 
parser.add_argument('-i', '--individuals',       type=int, help='Number of individuals') 
parser.add_argument('-r', '--temperarure_range', type=float, help='Range for testet temperatures (Temp +- (range/2))') 
parser.add_argument('--time_range',              type=int, help='Range for testet times        (Time +- (range/2))') 
parser.add_argument('--time_step',               type=int, help='Step between times') 
parser.add_argument('--strands',                 type=int, help='Number of strands') 
parser.add_argument('--strands_range',           type=int, help='Range of strands number') 
parser.add_argument('--strands_step',            type=int, help='step in strends number') 
parser.add_argument('-c', '--csv',               type=str, help='Store test results in csv file')
parser.add_argument('-p', '--plot',              type=str, help='Store ratio plot to file')
parser.add_argument('-s', '--samples',           type=int, help='How many samples per 1deg')
parser.add_argument('--target',                  type=str, help='Expected target output of simulation')

args = parser.parse_args()

target = args.target

if target is None:
    print("[INFO] No target specified using output from domain level simulation as target")
    target = subprocess.check_output([args.sdcsim, args.sdcasm, '-f', 'assembly'])
    target = target.decode("ascii").split("\n")
    target = target[-2]

if args.temperarure is None:
    print("[ERROR] No center temperature is defined using -T")
    exit(1)

if args.temperarure_range is None:
    print("[ERROR] No temperature range is definaded using -r")
    exit(1)

if args.individuals is None:
    print("[ERROR] number of individuals is not set using -i")
    exit(1)

if args.samples is None:
    print("[ERROR] number of samples per degree / second is not set using -s")
    exit(1)

if args.time is None:
    print("[ERROR] time of simulation must be definaded")
    exit(1)

if args.time_range is None:
    print("[ERROR] time range of simulation must be definaded")
    exit(1)

if args.time_step is None:
    print("[ERROR] time step of simulation must be definaded")
    exit(1)

if args.strands is None:
    print("[ERROR] number of strends in simulation must be definaded")
    exit(1)

if args.strands_range is None:
    print("[ERROR] range of number of strends in simulation must be definaded")
    exit(1)

if args.strands_step is None:
    print("[ERROR] step of number of strends in simulation must be definaded")
    exit(1)


print(f"[INFO] Going to run {args.sdcsim} {args.sdcasm} -n -f assembly -T {args.temperarure - args.temperarure_range/2}..{args.temperarure + args.temperarure_range/2}")

outputs = []

def forTemp(sdcsim, sdcasm, temp, individuals, target, time, strands):
    outputs = [temp, 0]
    for i in range(individuals):
        out = subprocess.check_output([sdcsim, sdcasm, '-n', '-f', 'assembly', '-T', str(temp), '-t', str(time), '-S', str(strands)])
        out = out.decode("ascii").split("\n")
        out = out[-2]
        if out == target:
            outputs[-1] += 1

    outputs[-1] /= individuals

    return outputs


fig, ax = plt.subplots()
Alloutputs = []

for strands in range(int(args.strands - args.strands_range / 2), int(args.strands + args.strands_range / 2), args.strands_step):
    for time in range(int(args.time - args.time_range / 2), int(args.time + args.time_range / 2), args.time_step):
        print(f"Time: {time} Strands: {strands}")
        outputs = Parallel(n_jobs=38)(delayed(forTemp)(args.sdcsim, args.sdcasm, ((i / args.samples - args.temperarure_range / 2) + args.temperarure), args.individuals, target, time, strands) for i in range(int(args.temperarure_range * args.samples)))

        if args.plot is not None:
            outputs = np.array(outputs)
            ax.plot(outputs[:,0], outputs[:,1], label=f"time: {time} strands: {strands}")

        for tempI in range(int(args.temperarure_range * args.samples)):
            Alloutputs.append([
                time,
                strands,
                outputs[tempI][0],
                outputs[tempI][1]
            ])


if args.plot is not None:

    ax.set(xlabel='Temperature C', ylabel='Target strands',
        title='Ratio of target strands')
    ax.grid()
    ax.legend()

    fig.savefig(args.plot)

if args.csv is not None:
    f = open(args.csv, "w+")
    f.write("time;strands;temp;finals\n")
    for output in Alloutputs:
        f.write(";".join(str(v) for v in output) + '\n')

    f.close()
