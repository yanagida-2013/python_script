#!/usr/bin/env python3
# coding: utf-8
#
# s.yanagida 2014
# 捕獲ガンマ線スペクトルにあわせて表記するレベルを
# gnuplotに適した形式で吐き出すプログラム。
# 元データとしてcconeのデータファイルを使っている。
#

import sys
import time

element = ['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', 'Na',
           'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', 'Sc', 'Ti',
           'V', 'Cr',  'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As',
           'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
           'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe', 'Cs',
           'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
           'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir',
           'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',
           'Ac', 'Th', 'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
           'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
           'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']


def read(z, a, neutron_energy):
    filename = ('levels/z%3d.lev' % z).replace(" ", "0")
    nucleus = str(a) + element[z-1]
    next_nucleus = str(a+1) + element[z-1]
    print(filename)
    print(nucleus)

    # read level files
    flag = False
    levels = []
    BnPlusEn = 0.0
    try:
        for line in open(filename, 'r'):
            if line[0:len(next_nucleus)] == next_nucleus:
                flag = False

            if flag and line[22] != ' ':
                # print(line.rsplit())
                levels.append(float(line.rsplit()[1]))

            if line[0:len(nucleus)] == nucleus:
                BnPlusEn = float(line.rsplit()[7]) + neutron_energy
                flag = True
    except IOError:
        print('File not found.')
        print('No levels are saved for %s' % nucleus)
        sys.exit()
    if len(levels) == 0:
        print('No levels were found for %s' % nucleus)
        sys.exit()
    else:
        print('%d levels found for %s' % (len(levels), nucleus))
        return levels, BnPlusEn, neutron_energy


def write(z, a, min, max, neutron_energy):
    levels = read(z, a, neutron_energy)
    filename = ('levels/z%3d.lev' % z).replace(" ", "0")
    nucleus = str(a) + element[z-1]
    outfile = nucleus + ".dat"
    out = open(outfile, 'w')
    out.write("# Levels of %s\n" % nucleus)
    out.write("# Generated for gnuplot input by barcode.py\n")
    out.write("# Rawdata file   :  %s\n" % filename)
    out.write("# Generated time :  %s\n" % time.ctime())
    out.write("# Neutron Energy :  %.3f[MeV]\n\n\n" % levels[2])

    for level in levels[0]:
        string = '%6.5e\t%6.5e\n%6.5e\t%6.5e\n\n\n' % (levels[1]-level, min, levels[1]-level, max)
        out.write(string)
    print('# In gnuplot script')
    print('plot "%s"w l' % outfile)
    out.close


def usage():
    print("""
Usage:
    ./barcode.py Z_num A_num y_min y_max Neutron_energy[MeV]
Example:
    ./barcode.py 56 139 10 100 0.03    #=> Ba139
    ./barcode.py 79 198 10 100 0.03    #=> Au198
    """)


def main():
    argvs = sys.argv
    if len(argvs) != 6:
        usage()
        sys.exit()

    try:
        z = int(argvs[1])
        if z < 1 or z > len(element):
            print("Z number is out of range")
            sys.exit
        a = int(argvs[2])
        min = float(argvs[3])
        max = float(argvs[4])
        neutron_energy = float(argvs[5])
        write(z, a, min, max, neutron_energy)
    except:
        print("Some inputs are wrong.")
        sys.exit()

if __name__ == "__main__":
    main()
