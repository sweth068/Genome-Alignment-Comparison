from prettytable import PrettyTable

mismatch_counts = {}

with open('/home/swetha/final/minimap_table_output/mini_negative_mismatch_distribution.txt', 'r') as file:
    for line in file:
        if line.startswith('Position'):
            _, mismatches_str = line.split(',')
            mismatches = int(mismatches_str.split(':')[1])

            if mismatches in mismatch_counts:
                mismatch_counts[mismatches] += 1
            else:
                mismatch_counts[mismatches] = 1

#Create a table
table = PrettyTable(['Mismatches', 'Count'])

#Add data to the table
for mismatches, count in mismatch_counts.items():
    table.add_row([mismatches, count])

#Print the table
print(table)
