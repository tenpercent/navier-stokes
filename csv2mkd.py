#!/usr/bin/env python3

results_csv_file = open('results.csv')
results_csv = results_csv_file.read()
results_csv_file.close()

headers = ('Время', '$\\omega$',
           '$\\tau$', '$h$',
           'невязка $V$ в $C$',
           'невязка $V$ в $L_2$',
           'невязка $G$ в $C$',
           'невязка $G$ в $L_2$')

widths = (8, 8, 11, 11, 17, 19, 17, 19)

def print_items(items):
    print(' | '.join((h + ' ' * (widths[i] - len(items[i])) for i, h in enumerate(items))))

print_items(headers)
print(' | '.join(('-' * width for width in widths)))

for line in results_csv.split('\n')[1:-1]:
    line = [el.strip() for el in line.split(',')]
    newline = ['%2.3f с' % float(line[1][:-2]),
               '%1.1f' % float(line[2])]
    newline += ['%.5e' % float(line[i]) for i in range(3, 9)]
    print_items(newline)
