#!/usr/bin/env python3

import sys
if len(sys.argv) > 1:
    results_csv_file = open(sys.argv[1])
    results_csv = results_csv_file.read()
    results_csv_file.close()
else:
    results_csv = sys.stdin.read()

headers = ('$\\tau$', '$h$',
           'Невязка $V$ в $C$',
           'Невязка $V$ в $L_2$',
           'Невязка $G$ в $C$',
           'Невязка $G$ в $L_2$',
           'Время')

widths = (6, 6, 17, 19, 17, 19, 7)

def print_items(items):
    print(' | '.join((h + ' ' * (widths[i] - len(items[i])) for i, h in enumerate(items))).strip())

print_items(headers)
print(' | '.join(('-' * width for width in widths)))

for line in results_csv.split('\n')[1:-1]:
    line = [el.strip() for el in line.split(',')]
    newline = [('%1.3f' % float(line[3])).rstrip('0'),
               ('%1.3f' % float(line[4])).rstrip('0')]
    newline += ['%.5e' % float(line[i]) for i in range(5, 9)]
    newline += ['%2.3f с' % float(line[1][:-2])]
    print_items(newline)
