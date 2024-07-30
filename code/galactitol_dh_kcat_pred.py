#@title Example prediction:
# prediction_type = "kcat" #@param ["kcat", "Km", "Ki"]
parameter = 'kcat' # @param ["kcat", "Km", "Ki"] {allow-input: true}
parameter = parameter.lower()

import time
import pandas as pd
import numpy as np
from IPython.display import Image, display
from rdkit import Chem
from rdkit.Chem.Draw.IPythonConsole import ShowMols

uniprot_id = "" #@param {type:"string"}
#@markdown If you do not have a uniprot-id, enter some name (for eg: "enzyme1")
sequence = 'MSFGSKAAERLANKIVLITGASSGIGAATAREFASAANGNIRLVLTARRESRLTELAKELSSTYAAIRVHTAKLDVSDTATIKPFISALPAEFADIDVLINNAGKALGKDVVGDISLDDISGMLQTNVLGLINVTQAVLPGMKARNSGSIVNIGSIAGREPYPGGSVYCASKAAVKSFSHALRKELISTRVRVLEVDPGAVETEFSVVRFHGDKAAADKVYEGTEPLTPEDIAEIIVFGVSRKENTVLAEALVFPSHQAGAVHVYKKSV' #@param {type:"string"}
SMILES = "NC(=O)C1=CN([C@H]2O[C@@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@@H](O)C3O)[C@@H](O)C2O)C=CC1.[H+].OC[C@@H](O)[C@H](O)[C@H](O)C(=O)CO.NC(=O)c1ccc[n+]([C@H]2O[C@@H](COP(=O)([O-])OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@@H](O)C3O)[C@@H](O)C2O)c1.OC[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)CO" #@param {type:"string"}

def create_csv_sh(parameter, uni, seq, smi):
  try:
    mol = Chem.MolFromSmiles(smi)
    smi = Chem.MolToSmiles(mol)
  except:
    print('Invalid SMILES input!')
    print('Correct your input! Exiting..')
    return
  valid_aas = list('ACDEFGHIKLMNPQRSTVWY')
  for aa in seq:
    if not aa in valid_aas:
      print('Invalid Enzyme sequence input!')
      print('Correct your input! Exiting..')
      return
  if parameter=='kcat':
    if '.' in smi:
      x = smi.split('.')
      y = sorted(x)
      smi = '.'.join(y)
  f = open(f'{uni}_{parameter}_input.csv', 'w')
  f.write('name,sequence,SMILES,pdbpath\n')
  f.write(f'{uni},{seq},{smi},{uni}.pdb\n')
  # f.write(f'{uni},{seq},{smi},{uni}.pdb')
  f.close()

  f = open(f'predict.sh', 'w')
  f.write(f'''
TEST_FILE_PREFIX={uni}_{parameter}
RECORDS_FILE=${{TEST_FILE_PREFIX}}.json
CHECKPOINT_DIR=./production_models/{parameter}/

python ./scripts/create_pdbrecords.py --data_file ${{TEST_FILE_PREFIX}}_input.csv --out_file ${{RECORDS_FILE}}
python predict.py --test_path ${{TEST_FILE_PREFIX}}_input.csv --preds_path ${{TEST_FILE_PREFIX}}_output.csv --checkpoint_dir $CHECKPOINT_DIR --uncertainty_method mve --smiles_column SMILES --individual_ensemble_predictions --protein_records_path $RECORDS_FILE
''')
  f.close()

  print('Input success!')
  print('Enzyme sequence length:', len(sequence))
  print('Substrate structure:')
  display(ShowMols([mol]))

  return seq, smi

seq, smi = create_csv_sh(parameter, uniprot_id, sequence, SMILES)

print('Predicting.. This will take a while..\n')

!export PROTEIN_EMBED_USE_CPU=1;./predict.sh >/dev/null 2>&1

from IPython.display import display, Latex, Math

def get_predictions(parameter, uniprot_id):
  confidence_limits = {'ki': np.array([0.6867484972083392,
                            0.9272913695545892,
                            1.3923101677384617,
                            5.328335398525296]),
                     'km': np.array([0.13853424987416604,
                            0.2213414848095112,
                            0.37158622530077795,
                            1.5143408432109369]),
                     'kcat': np.array([0.2920494204766375,
                              0.4801333784006298,
                              0.7838441054409631,
                              3.907686444864487])}
  confidence_limits = {key: np.power(10,val) for key, val in confidence_limits.items()}
  df = pd.read_csv(f'{uniprot_id}_{parameter}_output.csv')
  unit = ' mM'
  if parameter=='kcat':
    parameter_print = 'k_{cat}'
    parameter_print_log = 'log_{10}(k_{cat})'
    target_col = 'log10kcat_max'
    unit = ' s^{-1}'
  elif parameter=='km':
    target_col = 'log10km_mean'
    parameter_print = 'K_{m}'
    parameter_print_log = 'log_{10}(K_{m})'
  else:
    target_col = 'log10ki_mean'
    parameter_print = 'K_{i}'
    parameter_print_log = 'log_{10}(K_{i})'

  unc_col = f'{target_col}_mve_uncal_var'
  model_cols = [col for col in df.columns if col.startswith(target_col) and 'model_' in col]

  unc = df[unc_col].iloc[0]

  prediction = df[target_col].iloc[0]
  prediction_linear = np.power(10, prediction)

  model_out = df[target_col].iloc[0]
  model_outs = np.array([df[col].iloc[0] for col in model_cols])
  # print(model_outs)
  epi_unc = np.var(model_outs)#np.sum(np.power(2, model_outs))/10. - np.power(2, model_out)
  alea_unc = unc - epi_unc
  epi_unc = np.sqrt(epi_unc)
  alea_unc = np.sqrt(alea_unc)
  unc = np.sqrt(unc)
  ub = np.power(10, df[target_col].iloc[0]+3*unc)
  lb = np.power(10, df[target_col].iloc[0]-3*unc)

  # print(unc-epi_unc-alea_unc)
  # def display_outs(prediction_type, out, alea_output, epi_output, unit):
  print('Result:\n')
  display(Math((parameter_print + f' = {prediction_linear:.5f}'+ unit)))
  print('\n')
  display(Math((parameter_print_log + f' = {prediction:.5f}')))
  # print('UB: ', ub, 'LB: ', lb)
  limits = confidence_limits[parameter]
  # print(unc, limits[j])
  display(Math(('SD_{total}'+f' = {unc:.5f}')))
  display(Math(('SD_{aleatoric}'+f' = {alea_unc:.5f}')))
  display(Math(('SD_{epistemic}'+f' = {epi_unc:.5f}')))

get_predictions(parameter,uniprot_id)