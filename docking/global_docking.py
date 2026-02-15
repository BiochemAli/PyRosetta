#
#//  global_docking.py
#//  Pyrosetta_github
#//
#//  Created by A. Seltmann on 15.02.26.
#//
#//  adapted from     github.com/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/08.01-Ligand-Docking-XMLObjects.ipynb
#//

import logging
logging.basicConfig(level=logging.INFO)
import matplotlib
import os
import pandas as pd
import pyrosetta
import pyrosetta.distributed.viewer as viewer
import seaborn
seaborn.set()
import sys
import argparse

parser = argparse.ArgumentParser(description="PyRosetta global docking w/ XML")

parser.add_argument('-i', '--input_pdb', required=True, help='Path to PDB file')
parser.add_argument('-p', '--params_file', required=True, help='Path to parameter file of the ligand')
parser.add_argument('-n', '--number_outputs', required=True, type=int, default=10, help='Number of generated models')
args = parser.parse_args()

ligand_params = args.params_file
flags = f"""
-ignore_unrecognized_res 1
-extra_res_fa {ligand_params}
-ex1 -ex2
-use_input_sc
-no_optH false
-flip_HNQ
"""
pyrosetta.distributed.init(flags)
pose = pyrosetta.io.pose_from_file(filename=args.input_pdb)
scorefxn = pyrosetta.create_score_function("ref2015")

xml = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="fa_standard" weights="ref2015.wts"/>
  </SCOREFXNS>
  <RESIDUE_SELECTORS>
    <Chain name="chX" chains="X"/>
  </RESIDUE_SELECTORS>
  <SIMPLE_METRICS>
    <RMSDMetric name="rmsd_chX" residue_selector="chX" reference_name="store_native" residue_selector_ref="chX" robust="true" rmsd_type="rmsd_all" />
  </SIMPLE_METRICS>
  <SCORINGGRIDS ligand_chain="X" width="25">
    <ClassicGrid grid_name="vdw" weight="1.0"/>
  </SCORINGGRIDS>
  <LIGAND_AREAS>
    <LigandArea name="docking_sidechain_X" chain="X" cutoff="6.0" add_nbr_radius="true" all_atom_mode="true" minimize_ligand="10"/>
    <LigandArea name="final_sidechain_X" chain="X" cutoff="6.0" add_nbr_radius="true" all_atom_mode="true"/>
    <LigandArea name="final_backbone_X" chain="X" cutoff="7.0" add_nbr_radius="false" all_atom_mode="true" Calpha_restraints="0.3"/>
  </LIGAND_AREAS>
  <INTERFACE_BUILDERS>
    <InterfaceBuilder name="side_chain_for_docking" ligand_areas="docking_sidechain_X"/>
    <InterfaceBuilder name="side_chain_for_final" ligand_areas="final_sidechain_X"/>
    <InterfaceBuilder name="backbone" ligand_areas="final_backbone_X" extension_window="3"/>
  </INTERFACE_BUILDERS>
  <MOVEMAP_BUILDERS>
    <MoveMapBuilder name="docking" sc_interface="side_chain_for_docking" minimize_water="true"/>
    <MoveMapBuilder name="final" sc_interface="side_chain_for_final" bb_interface="backbone" minimize_water="true"/>
  </MOVEMAP_BUILDERS>
  <MOVERS>
    <SavePoseMover name="spm" restore_pose="0" reference_name="store_native"/>
    <Transform name="transform" chain="X" box_size="20.0" move_distance="10" angle="360" initial_perturb="2" cycles="500" repeats="5" temperature="1000"/>
    <HighResDocker name="high_res_docker" cycles="9" repack_every_Nth="3" scorefxn="fa_standard" movemap_builder="docking"/>
    <FinalMinimizer name="final" scorefxn="fa_standard" movemap_builder="final"/>
  </MOVERS>
  <FILTERS>
      <LigInterfaceEnergy name="interfE" scorefxn="fa_standard" energy_cutoff="0.0" confidence="0"/>
      <SimpleMetricFilter name="rmsd_chX" metric="rmsd_chX" cutoff="999999." comparison_type="lt" confidence="0"/>
  </FILTERS>
  <PROTOCOLS>
    <Add mover="spm"/>
    <Add mover="transform"/>
    <Add mover="high_res_docker"/>
    <Add mover="final"/>
    <Add filter="interfE"/>
    <Add filter="rmsd_chX"/>
  </PROTOCOLS>
</ROSETTASCRIPTS>
""").get_mover("ParsedProtocol")

if not os.getenv("DEBUG"):
    working_dir = os.getcwd()
    output_dir = "outputs"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    os.chdir(output_dir)

    jd = pyrosetta.toolbox.py_jobdistributor.PyJobDistributor(pdb_name=args.input_pdb,
                                                              nstruct=args.number_outputs,
                                                              scorefxn=scorefxn)
    jd.native_pose = pose
    
    df_list = []
    
    while not jd.job_complete:
        test_pose = pose.clone()
        xml.apply(test_pose)
        test_df = pd.DataFrame.from_records(dict(test_pose.scores), index=[jd.current_name])
        df_list.append(test_df)
        jd.output_decoy(test_pose)
    
    df = pd.concat(df_list, ignore_index=False)
    
    os.chdir(working_dir)
