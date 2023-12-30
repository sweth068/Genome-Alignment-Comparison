mismatches = {}

with open('mini_negative_pileup.txt', 'r') as pileup_file:
    for line in pileup_file:
        fields = line.strip().split('\t')
        position = int(fields[1])
        bases = fields[4]
        qualities = fields[5]

        #Count mismatches (base not matching the reference)
        mismatch_count = sum(1 for base, qual in zip(bases, qualities) if base.upper() not in 'ACGT' or ord(qual) - 33 < 20)

        mismatches[position] = mismatch_count

#Write mismatch distribution to a new file
output_file = '/home/swetha/final/minimap_table_output/mini_negative_mismatch_distribution.txt'
with open(output_file, 'w') as output:
    for position, count in sorted(mismatches.items()):
        output.write(f'Position: {position}, Mismatches: {count}\n')

print(f'Mismatch distribution written to {output_file}')
