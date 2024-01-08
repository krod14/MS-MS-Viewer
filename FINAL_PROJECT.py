# import modules and functions

import sys
import gzip
from base64 import b64decode
from array import array
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

# file command-line argument and namespace

try:
    file = sys.argv[1]
except IndexError:
    print("Please input an mzxml.gz file.")
    sys.exit()

ns = '{http://sashimi.sourceforge.net/schema/}'

# scan and peptide command-line arguments

try:
    scan_number = sys.argv[2]
except IndexError:
    print("Please input a scan number.")
    sys.exit()

try:
    peptide = sys.argv[3]
    peptide = peptide.upper()
except IndexError:
    print("Please input a peptide sequence.")
    sys.exit()

# calculate molecular weight

mw = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
      'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
      'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
      'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }

def MolWt(pep):
    molwt = 0
    for aa in pep:
        aawt = mw.get(aa)
        molwt += aawt
    return molwt

# calculate m/z b-ion, call function

def mzb_ion(pep):
    b_ion = []
    mz = []
    for aa in pep:
        b_ion.append(aa)
        result = MolWt(b_ion) + 1
        mz.append(result)
    mz_b_d = dict((f'b{index}', value) for index, value in enumerate(mz, start=1))
    return mz_b_d

try:
    b_ion = mzb_ion(peptide)
except TypeError:
    print("Please input a valid peptide sequence.")
    sys.exit()

# calculate m/z y-ion, call function

def mzy_ion(pep):
    pep = ''.join(reversed(pep))
    y_ion = []
    mz = []
    for aa in pep:
        y_ion.append(aa)
        result = MolWt(y_ion) + 19
        mz.append(result)
    mz_y_d = dict((f'y{index}', value) for index, value in enumerate(mz, start=1))
    return mz_y_d

y_ion = mzy_ion(peptide)

# access peaks element, put (mzs, ints) into list as tuples

try:
    with gzip.open(file) as gzipped_file:
        for event,ele in ET.iterparse(gzipped_file):
            if ele.tag == ns+'scan':
                if ele.attrib['num'] == scan_number:
                    peakselt = ele.find(ns+'peaks')
                    peaks = array('f',b64decode(peakselt.text))
                    if sys.byteorder != 'big':
                       peaks.byteswap()
                    mzs = peaks[::2]
                    ints = peaks[1::2]
except FileNotFoundError:
    print("Please input a valid mzxml.gz file.")
    sys.exit()
    
try:    
    mzs_ints_l = list(zip(mzs, ints))
except NameError:
    print("Please input a valid scan number.")
    sys.exit()

# find max intensity, change intensity in tuple to relative absorbance

max_int = max(ints[1] for ints in mzs_ints_l)

mzs_relabs_l = [(mzs, round((ints / max_int) * 100, 3)) for mzs, ints in mzs_ints_l]

# matching peaks? if so, add to dictionary with matching ion as key and (mz, int) as tuple value

b_y_match_d = {}

for mzs, relabs in mzs_relabs_l:
    for i, mz in b_ion.items():
        lower_b = mz - 0.3
        upper_b = mz + 0.3
        if lower_b <= mzs <=  upper_b:
            b_y_match_d[i] = (mzs,relabs)
    for i, mz in y_ion.items():
        lower_y = abs(mz) - 0.3
        upper_y = abs(mz) + 0.3
        if lower_y <= mzs <=  upper_y:
            b_y_match_d[i] = (mzs,relabs)

# pick apart information, annotate, and plot

all_mzs, all_relabs = zip(*mzs_relabs_l)
ion_keys, mz_relabs_values = zip(*b_y_match_d.items())
match_mzs, match_relabs = zip(*mz_relabs_values)

match_mzs_d = {}

for ion, match_mzs, match_relabs in zip(ion_keys, match_mzs, match_relabs):
    if ion.startswith('b'):
        label_color = "blue"
    else:
        label_color = "red"
    if match_mzs in match_mzs_d:
        adjust = 2.5 + (match_mzs_d[match_mzs] - 1) * 2
        match_mzs_d[match_mzs] += 1
    else:
        adjust = 0
        match_mzs_d[match_mzs] = 1
    plt.text(match_mzs, match_relabs + adjust, ion, fontsize = 12, color = label_color, ha = 'center', va = 'bottom')

plt.stem(all_mzs, all_relabs, linefmt = 'k-', markerfmt = ' ', use_line_collection = True)
plt.xlabel('m/z', fontsize = 12)
plt.ylabel('Relative Abundance (%)', fontsize = 12)
plt.ylim(bottom = 0)
plt.suptitle('m/z vs. Relative Abundance (%)', x = 0.52, y = 0.97, fontsize = 15)
plt.title(scan_number + '  ' + peptide, fontsize = 10)

plt.show()

