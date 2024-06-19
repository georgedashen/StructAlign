import os, sys
import glob

infolder = sys.argv[1]
outfile = sys.argv[2]

# The suffix of the files you want to find
suffix = '.tmscore'  # Replace with your specific suffix

# Find all files with the specified suffix, including in subdirectories
files_to_concatenate = glob.glob(f'{infolder}/**/*tmscore', recursive=True) #eg. Malisam, Malidup

# Open the output file in write mode
with open(outfile, 'w') as f:
    # Iterate over each file to concatenate
    for i, filename in enumerate(sorted(files_to_concatenate)):
        # Open each file in read mode
        with open(filename, 'r') as infile:
            # Read the content and write it to the output file
            if i == 0:
                f.write(infile.read())
            else:
                f.write('\n'.join(infile.read().splitlines()[1:])+'\n')

