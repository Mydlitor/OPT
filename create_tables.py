import csv

# Read the three test CSV files
test_files = ['lab5_test_a1.csv', 'lab5_test_a10.csv', 'lab5_test_a100.csv']
a_values = [1, 10, 100]

# Create Table 1 - Test functions (all a values combined)
with open('lab5_table1.csv', 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    
    # Write header
    writer.writerow(['a', 'w', 'x1_start', 'x2_start', 'x1_opt', 'x2_opt', 'f1_opt', 'f2_opt', 'f_opt', 'f_calls'])
    
    # Process each test file
    for a, filename in zip(a_values, test_files):
        with open(filename, 'r') as infile:
            reader = csv.DictReader(infile)
            for row in reader:
                writer.writerow([
                    a,
                    row['w'],
                    row['x1_start'],
                    row['x2_start'],
                    row['x1_opt'],
                    row['x2_opt'],
                    row['f1_opt'],
                    row['f2_opt'],
                    row['f_opt'],
                    row['f_calls']
                ])

print("Created lab5_table1.csv with combined test function results")

# Create Table 2 - Beam problem (just copy and rename)
import shutil
shutil.copy('lab5_beam.csv', 'lab5_table2.csv')
print("Created lab5_table2.csv from beam results")

# Show summary
print(f"\nTable 1 rows: {sum(1 for line in open('lab5_table1.csv')) - 1}")
print(f"Table 2 rows: {sum(1 for line in open('lab5_table2.csv')) - 1}")
