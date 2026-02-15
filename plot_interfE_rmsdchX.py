#
#//  global_docking.py
#//  Pyrosetta_github
#//
#//  Created by A. Seltmann on 15.02.26.
#//
#//  adapted from     github.com/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/08.01-Ligand-Docking-XMLObjects.ipynb
#//

import glob
import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import argparse

parser = argparse.ArgumentParser(description="PyRosetta plotting interfE against rmsd_chX")

parser.add_argument('-p', '--params_file', required=True, help='Path to parameter file of the ligand')
args = parser.parse_args()

df = pd.read_json('outputs/*.fasc', lines=True)
df_sorted = df.sort_values(by='interfE')

# Für neuere Pandas-Versionen (append ist deprecated)
if not os.getenv("DEBUG"):
    pdb_filenames = glob.glob("outputs/*.pdb")
    ligand_params = args.params_file
    native_pdb_filename = f"outputs/{df_sorted['filename'].head(1)}"

    flags = f"""
    -extra_res_fa {ligand_params} 
    -in:file:native {native_pdb_filename}
    -ignore_unrecognized_res 1 
    -mute all
    """
    pyrosetta.distributed.init(flags)
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
        <RMSDMetric name="rmsd_chX" use_native="true" residue_selector="chX" residue_selector_ref="chX" robust="true" rmsd_type="rmsd_all" />
      </SIMPLE_METRICS>
      <FILTERS>
          <LigInterfaceEnergy name="interfE" scorefxn="fa_standard" energy_cutoff="0.0" confidence="0"/>
          <SimpleMetricFilter name="rmsd_chX" metric="rmsd_chX" cutoff="999999." comparison_type="lt" confidence="0"/>
      </FILTERS>
      <PROTOCOLS>
        <Add filter="interfE"/>
        <Add filter="rmsd_chX"/>
      </PROTOCOLS>
    </ROSETTASCRIPTS>
    """).get_mover("ParsedProtocol")

    df_list = []
    for pdb_filename in pdb_filenames:
        test_pose = pyrosetta.io.pose_from_file(filename=pdb_filename)
        xml.apply(test_pose)
        test_df = pd.DataFrame.from_records(dict(test_pose.scores), index=[pdb_filename.split("/")[-1]])
        df_list.append(test_df)
    
    df = pd.concat(df_list, ignore_index=False)

    matplotlib.rcParams['figure.figsize'] = [6.0, 4.0]
    seaborn.scatterplot(x="rmsd_chX", y="interfE", data=df)
    plt.show()

