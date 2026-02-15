#
#//  fastrelax.py
#//  Pyrosetta_github
#//
#//  Created by A. Seltmann on 15.02.26.
#//

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

parser = argparse.ArgumentParser(description='PyRosetta FastRelax')

parser.add_argument('-i', '--input_pdb', required=True, help='Path to PDB file')
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

relax.apply(testPose)

testPose.dump_pdb(f"{args.input_pdb[:-4]}_relaxed.pdb")
print(f"Done! Saved as {args.input_pdb[:-4]}_relaxed.pdb")
