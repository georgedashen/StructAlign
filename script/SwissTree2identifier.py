import toytree
import os

os.chdir('examples/swiss_trees_updated')
for t in os.listdir('.'):
    if not t.endswith('nhx'):
        continue
    tree = toytree.tree(t)
    nodedict = tree.get_node_dict()
    entries = [i for i in list(nodedict.values()) if i not in ['None','failed']]
    folder = t.split('_')[0]
    os.system(f'mkdir -p {folder}')
    os.system(f'cp {t} {folder}/{t}')
    with open(f'../../../SwissTree/{folder}/identifiers.txt','w') as f:
        for entry in entries:
            f.write(f'{entry}\n')
