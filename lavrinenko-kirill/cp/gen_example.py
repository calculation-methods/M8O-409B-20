import argparse
from pprint import pprint

from numpy import arange


def gen_example():
    map = dict()

    afind = []
    idx = 0
    for i in arange(0.1, 1, 0.1):
        for j in arange(0.1, 0.6, 0.1):
            if (round(i, 1), round(j, 1)) not in map:
                map[(round(i, 1), round(j, 1))] = idx
                idx += 1
            afind.append([round(i, 1), round(j, 1)])
    l = len(afind)

    # pprint(map)

    dx = 0.1
    table = [[0 for _ in range(l)] for _ in range(l)]
    val = [0 for _ in range(l)]
    for i, x in enumerate(afind):

        table[i][map[(x[0], x[1])]] = 4
        if x[0] - dx != 0:
            table[i][map[(round(x[0] - dx, 1), x[1])]] = -1
        else:
            val[i] = 100
        if x[1] - dx != 0:
            table[i][map[(x[0], round(x[1] - dx, 1))]] = -1
        if x[0] + dx != 1:
            table[i][map[(round(x[0] + dx, 1), x[1])]] = -1
        if x[1] + dx != 1:
            if round(x[1] + dx, 1) > 0.5:
                table[i][map[(x[0], 1 - round(x[1] + dx, 1))]] += -1
            else:
                table[i][map[(x[0], round(x[1] + dx, 1))]] = -1

        print(f"$$4T_{{{int(x[0]*10)},{int(x[1]*10)}}} - " +
              (f"T_{{{int(round(x[0] - dx, 1)*10)},{int(x[1]*10)}}}" if x[0] - dx != 0 else f"100") + " - " +
              (f"T_{{{int(x[0]*10)},{int(round(x[1] - dx, 1)*10)}}}" if x[1] - dx != 0 else f"0") + " - " +
              (f"T_{{{int(round(x[0] + dx, 1)*10)},{int(x[1]*10)}}}" if x[0] + dx != 1 else f"0") + " - " +
              (f"T_{{{int(x[0]*10)},{int(round(x[1] + dx, 1)*10)}}}" if x[1] + dx != 1 else f"0") + " = 0$$")

    return table, val


def write_to(table, val, filename):
    with open(filename, 'w') as f:
        f.write(f'{len(table)}\n')
        for row in table:
            for x in row:
                f.write(f'{x} ')
            f.write(f'\n')
        for x in val:
            f.write(f'{x} ')
        f.write(f'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, help='Output file')
    args = parser.parse_args()

    table, val = gen_example()
    write_to(table, val, args.output)


if __name__ == '__main__':
    main()
