import os,sys
import itertools
import pandas as pd

list1 = glob.glob('structure_data/pdb70_and_scope/*pdb')
list2 = glob.glob('structure_data/scope140_pdb/*ent')

list1 = [os.path.basename(i) for i in list1]
list2 = [os.path.basename(i) for i in list2]

# Generate all combinations
combinations = itertools.product(list1, list2)

# Create a DataFrame from the combinations
df = pd.DataFrame(combinations, columns=['Query', 'Target'])

df.to_csv('comparison.list', index=False)
