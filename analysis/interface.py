#!/usr/bin/env python3

import argparse
import os
from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

init('''
     -ex1 -ex2
     -use_input_sc
     -no_optH false
     -flip_HNQ
     -ignore_unrecognized_res
     -mute all
     '''
     )

parser = argparse.ArgumentParser(description='PyRosetta Interface Analysis w/ optional FastRelax')

parser.add_argument('-i', '--input_pdb', required=True, help='Path to PDB file')
parser.add_argument('-c', '--chains', required=True, help='Interface chains (e.g. A_B)')
parser.add_argument('-r','--relax', action='store_true', help='FastRelax (default: false)')
args = parser.parse_args()

pose = pose_from_pdb(args.input_pdb)
testPose = Pose()
testPose.assign(pose)
print("Pose loaded")

mm = MoveMap()
mm.set_bb(True)
mm.set_chi(True)

relax = FastRelax()
scorefxn = get_fa_scorefxn()
relax.set_scorefxn(scorefxn)
relax.set_movemap(mm)
relax.constrain_relax_to_start_coords(True)

if args.relax:
    print("Start FastRelax")
    relax.apply(testPose)
    print("FastRelax finished")

interface = InterfaceAnalyzerMover(args.chains)
interface.set_scorefunction(scorefxn)
print("Analyze interface.")
interface.apply(testPose)

if args.relax:
    testPose.dump_pdb(f"{args.input_pdb[:-4]}_relaxed_analyzed.pdb")
    print(f"Done! Saved as {args.input_pdb[:-4]}_relaxed_analyzed.pdb")
    with open(f'{args.input_pdb[:-4]}_relaxed_analyzed.pdb', 'r') as f:
        lines = f.readlines()
        print(''.join(lines[-21:]))
else:
    testPose.dump_pdb(f"{args.input_pdb[:-4]}_analyzed.pdb")
    print(f"Done! Saved as {args.input_pdb[:-4]}_analyzed.pdb")
    with open(f'{args.input_pdb[:-4]}_analyzed.pdb', 'r') as f:
        lines = f.readlines()
        print(''.join(lines[-21:]))
    
