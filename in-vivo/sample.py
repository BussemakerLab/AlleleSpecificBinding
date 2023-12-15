import sys

import pyprobound

# Read in file
_, input_file, bound_file, output_file = sys.argv
df = pyprobound.get_dataframe([input_file, bound_file])
df = df[df.index.str.fullmatch("[ACGT]*")]

# Get middle million lengths if pair-ended
lengths = df.index.str.len()
sort_lengths = lengths.sort_values()
if sort_lengths[0] != sort_lengths[-1]:
    min_len = sort_lengths[(len(sort_lengths) // 2) - 500_000]
    max_len = sort_lengths[(len(sort_lengths) // 2) + 500_000]
    df = df[(lengths >= min_len) & (lengths <= max_len)]

# Sample rows
df = df.sample(n=1_000_000, random_state=0)

# Write to file
df.to_csv(output_file, header=False, sep="\t")
