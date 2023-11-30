#!/usr/bin/env python

#  This is a sample script, with hard-wired variable names.  Copy this file
#   to your working directory, then modify the RDB file and column names,
#   and adjust the xlabal and ylabel names to match your project

#   Note that your PYTHONPATH variable must input parmed (most likely from
#   AmberTools), and your python must be able to import matplotlib, 
#   numpy, pandas and scipy.

from parmed.residue import AminoAcidResidue 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import pandas as pd
from scipy.stats import linregress

# This is how you would parse a RDB file. The data format line needs to be
# skipped by pandas, and skiprows=[1] means skip the second line 
# (first line is 0).
# If comments are present, skiprows needs to be adjusted to specify the correct
# absolute line

tab = pd.read_table('1d3z1.rdb', header=0, skiprows=[3], comment='#')

# This reads in another RDB and merges it with the above on the residue number,
# residue name, and atom name.
#    (note: the residue name is not needed for merging, but somehow that info
#     gets lost if not included here)

tab = tab.merge(pd.read_table('bmr17769.rdb', header=0, skiprows=[2],
                              comment='#'),
                on=['res', 'resname', 'atomname'])

# These are the names of the columns in the header of the RDB
CALC_HEADER = '1d3z1s'
EXPT_HEADER = 'exp'

# Name of the image file to write. Format is determined by filename extension
OUTPUT_FIGURE_NAME = '1d3z1s.pdf'

# These are the types of protons and carbons we want to plot
PROTONS = ['HA', '-CH-', '-CH2-', '-CH3', 'Aromatic-CH', 'HN',
           '-NH2', 'Aromatic-NH']
CARBONS = ['CA', 'Amide-C', '-HC-', '-H2C-', 'H3C-', 'Aromatic-C']

# This is the atom type map that maps the atom/residue pair names to the type of
# nucleus it is based on chemical environment
PROTON_TYPES = {
        '-CH3' : [('ALA', 'HB1'), ('ALA', 'HB2'), ('ALA', 'HB3'),
                  ('ILE', 'HG21'), ('ILE', 'HG22'), ('ILE', 'HG23'),
                  ('ILE', 'HD11'), ('ILE', 'HD12'), ('ILE', 'HD13'),
                  ('LEU', 'HD11'), ('LEU', 'HD12'), ('LEU', 'HD13'),
                  ('LEU', 'HD21'), ('LEU', 'HD22'), ('LEU', 'HD23'),
                  ('MET', 'HE1'), ('MET', 'HE2'), ('MET', 'HE3'),
                  ('THR', 'HG21'), ('THR', 'HG22'), ('THR', 'HG23'),
                  ('VAL', 'HG11'), ('VAL', 'HG12'), ('VAL', 'HG13'),
                  ('VAL', 'HG21'), ('VAL', 'HG22'), ('VAL', 'HG23')],

        '-CH2-' : [('ARG', 'HB2'), ('ARG', 'HB3'), ('ARG', 'HG2'),
                   ('ARG', 'HG3'), ('ARG', 'HD2'), ('ARG', 'HD3'),
                   ('ASN', 'HB2'), ('ASN', 'HB3'), ('ASP', 'HB2'),
                   ('ASP', 'HB3'), ('CYS', 'HB2'), ('CYS', 'HB3'),
                   ('GLN', 'HB2'), ('GLN', 'HB3'), ('GLN', 'HG2'),
                   ('GLN', 'HG3'), ('GLU', 'HB2'), ('GLU', 'HB3'),
                   ('GLU', 'HG2'), ('GLU', 'HG3'), ('HIS', 'HB2'),
                   ('HIS', 'HB3'), ('ILE', 'HG12'), ('ILE', 'HG13'),
                   ('LEU', 'HB2'), ('LEU', 'HB3'), ('LYS', 'HB2'),
                   ('LYS', 'HB3'), ('LYS', 'HG2'), ('LYS', 'HG3'),
                   ('LYS', 'HD2'), ('LYS', 'HD3'), ('LYS', 'HE2'),
                   ('LYS', 'HE3'), ('MET', 'HB2'), ('MET', 'HB3'),
                   ('MET', 'HG2'), ('MET', 'HG3'), ('PHE', 'HB2'),
                   ('PHE', 'HB3'), ('PRO', 'HB2'), ('PRO', 'HB3'),
                   ('PRO', 'HG2'), ('PRO', 'HG3'), ('PRO', 'HD2'),
                   ('PRO', 'HD3'), ('SER', 'HB2'), ('SER', 'HB3'),
                   ('TRP', 'HB2'), ('TRP', 'HB3'), ('TYR', 'HB2'),
                   ('TYR', 'HB3')],

        '-CH-' : [('ILE', 'HB'), ('LEU', 'HG'), ('THR', 'HB'), ('VAL', 'HB')],

        '-NH2' : [('ARG', 'HH21'), ('ARG', 'HH22'), ('ARG', 'HH11'),
                  ('ARG', 'HH12'), ('ASN', 'HD21'), ('ASN', 'HD22'),
                  ('GLN', 'HE21'), ('GLN', 'HE22'), ('LYS', 'HZ1'),
                  ('LYS', 'HZ2'), ('LYS', 'HZ3'),],

        '-SH/-OH' : [('CYS', 'HG'), ('SER', 'HG'), ('THR', 'HG1'),
                     ('TYR', 'HH')],

        'Aromatic-CH' : [('HIS', 'HD2'), ('HIS', 'HE1'), ('PHE', 'HD1'),
                         ('PHE', 'HD2'), ('PHE', 'HE1'), ('PHE', 'HE2'),
                         ('PHE', 'HZ'), ('TRP', 'HD1'), ('TRP', 'HE3'),
                         ('TRP', 'HZ3'), ('TRP', 'HH2'), ('TRP', 'HZ2'),
                         ('TYR', 'HD1'), ('TYR', 'HD2'), ('TYR', 'HE1'),
                         ('TYR', 'HE2')],

        'Aromatic-NH' : [('HIS', 'HE2'), ('HIS', 'HD1'), ('TRP', 'HE1')],
}

CARBON_TYPES = {
        'H3C-' : [('ALA', 'CB'), ('ILE', 'CD1'), ('ILE', 'CG2'),
                  ('LEU', 'CD1'), ('LEU', 'CD2'), ('MET', 'CE'),
                  ('THR', 'CG2'), ('VAL', 'CG1'), ('VAL', 'CG2')],

        '-H2C-' : [('ARG', 'CB'), ('ARG', 'CG'), ('ARG', 'CD'), ('ASN', 'CB'),
                   ('ASP', 'CB'), ('GLN', 'CB'), ('GLN', 'CG'), ('GLU', 'CB'),
                   ('GLU', 'CG'), ('HIS', 'CB'), ('ILE', 'CG1'), ('LYS', 'CB'),
                   ('LYS', 'CG'), ('LYS', 'CD'), ('LYS', 'CE'), ('MET', 'CB'),
                   ('MET', 'CG'), ('PHE', 'CB'), ('PRO', 'CB'), ('PRO', 'CG'),
                   ('PRO', 'CD'), ('SER', 'CB'), ('TRP', 'CB'), ('TYR', 'CB'),
                   ('LEU', 'CB')],

        '-HC-' : [('ILE', 'CB'), ('LEU', 'CG'), ('THR', 'CB'), ('VAL', 'CB')],

        'C=O' : [('ASN', 'CG'), ('ASP', 'CG'), ('GLN', 'CD'), ('GLU', 'CD')],

        'Aromatic-C' : [('HIS', 'CG'), ('HIS', 'CD2'), ('HIS', 'CE1'),
                        ('PHE', 'CG'), ('PHE', 'CD1'), ('PHE', 'CD2'),
                        ('PHE', 'CE1'), ('PHE', 'CE2'), ('PHE', 'CZ'),
                        ('TRP', 'CG'), ('TRP', 'CD1'), ('TRP', 'CE2'),
                        ('TRP', 'CD2'), ('TRP', 'CE3'), ('TRP', 'CZ3'),
                        ('TRP', 'CH2'), ('TRP', 'CZ2'), ('TYR', 'CG'),
                        ('TYR', 'CD1'), ('TYR', 'CD2'), ('TYR', 'CE2'),
                        ('TYR', 'CZ'), ('TYR', 'CE1')]
}

# Now build the CA, HA, N, and H entries from all of the residues

PROTON_TYPES['HA'] = list()
PROTON_TYPES['HN'] = list()
CARBON_TYPES['CA'] = list()
CARBON_TYPES['Amide-C'] = list()

for res in AminoAcidResidue.all_residues:
    PROTON_TYPES['HA'].append((res.abbr, 'HA'))
    PROTON_TYPES['HN'].append((res.abbr, 'H'))
    CARBON_TYPES['CA'].append((res.abbr, 'CA'))
    CARBON_TYPES['Amide-C'].append((res.abbr, 'C'))
PROTON_TYPES['HA'].append(('GLY', 'HA2'))
PROTON_TYPES['HA'].append(('GLY', 'HA3'))

# Below here controls the appearance of the plots. Above is where most of the
# modifications will be needed
legendfont = FontProperties()
legendfont.set_size(20)

symbols = ['o', 'x', '^', '*', 'h', '+', 'D', '8', 'v', '<', '>']
LS = len(symbols)
colors = ['r', 'b', 'g', 'k', 'm']
LC = len(colors)

tab['diff'] = tab[CALC_HEADER] - tab[EXPT_HEADER]

alltypes = dict()
alltypes.update(PROTON_TYPES)
alltypes.update(CARBON_TYPES)

def map_to_type(series):
    key = (series.resname, series.atomname)
    for type, items in alltypes.items():
        if key in items: return type
    if series.atomname == 'N':
        return 'N'
    return 'Unknown'

tab['chemtype'] = tab.apply(map_to_type, axis=1)

fig = plt.figure(0, figsize=(15,15))

ax = fig.add_subplot(111)
ax.tick_params(labelsize=23)
for key in ax.spines:
    ax.spines[key].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off',
               labelsize=23)
ax.set_xlabel('Experimental chemical shifts (ppm)', fontdict=dict(fontsize=30))
ax.set_ylabel('AF-QM/MM calculated chemical shifts (ppm)',
              fontdict=dict(fontsize=30))

# Add the hydrogens
ax = fig.add_subplot(221)
ax.tick_params(labelsize=23)
pls = []
names = []
for i, proton in enumerate(PROTONS):
    if proton in ('HN', '-NH2', 'Aromatic-NH'): continue
    marker = symbols[i%LS]
    color = colors[i%LC]
    data = tab.loc[tab.chemtype == proton]

    pl, = ax.plot(data[EXPT_HEADER], data[CALC_HEADER], marker=marker,
                  color=color, markersize=8, ls='None')
    pls.append(pl)
    names.append(proton)

def findh(s):
    return s.chemtype in PROTONS and \
            s.chemtype not in ['HN', '-NH2', 'Aromatic-NH']

all_hs = tab.loc[tab.apply(findh, axis=1)]
m, b, r, p, err = linregress(all_hs[EXPT_HEADER], all_hs[CALC_HEADER])
y = lambda x: m*x+b
hmin = min(all_hs[EXPT_HEADER].min(), all_hs[CALC_HEADER].min())
hmax = max(all_hs[EXPT_HEADER].max(), all_hs[CALC_HEADER].max())
mse = ((all_hs[CALC_HEADER] - all_hs[EXPT_HEADER])**2).mean()
ax.text(0.34, 0.08, '$R = %.3f$\n$RMSE = %.2f\ ppm$' % (r, np.sqrt(mse)),
        fontsize=30, transform=ax.transAxes)

ax.plot([-2, 10], [-2, 10], 'k-', linewidth=2)
ax.plot([hmin, hmax], [y(hmin), y(hmax)], 'k--', linewidth=2)
ax.legend(pls, names, numpoints=1, loc='upper left', prop=legendfont)
ax.set_title('$^1$H Shifts', fontdict=dict(fontsize=30))

# Add the alpha-carbons
ax = fig.add_subplot(222)
ax.tick_params(labelsize=23)
ax.set_title(r'$^{13}$C$\alpha$ Shifts', fontdict=dict(fontsize=30))
data = tab.loc[tab.atomname == 'CA']
ax.plot(data[EXPT_HEADER], data[CALC_HEADER], marker='o', color='k',
        markersize=8, ls='None')
ax.plot([40, 70], [40, 70], 'k-', linewidth=2)

m, b, r, p, err = linregress(data[EXPT_HEADER], data[CALC_HEADER])
y = lambda x: m*x+b
hmin = min(data[EXPT_HEADER].min(), data[CALC_HEADER].min())
hmax = max(data[EXPT_HEADER].max(), data[CALC_HEADER].max())
mse = ((data[EXPT_HEADER] - data[CALC_HEADER])**2).mean()
ax.text(0.34, 0.08, '$R = %.3f$\n$RMSE = %.2f\ ppm$' % (r, np.sqrt(mse)),
        fontsize=30, transform=ax.transAxes)
ax.plot([hmin, hmax], [y(hmin), y(hmax)], 'k--', linewidth=2)

# Add the rest of the carbons
ax = fig.add_subplot(223)
ax.tick_params(labelsize=23)
pls = []
names = []
for i, carbon in enumerate(CARBONS):
    if carbon in ('CA'): continue
    marker = symbols[i%LS]
    color = colors[i%LC]
    data = tab.loc[tab.chemtype == carbon]

    pl, = ax.plot(data[EXPT_HEADER], data[CALC_HEADER], marker=marker,
                  color=color, markersize=8, ls='None')
    pls.append(pl)
    names.append(carbon)

def findc(s):
    return s.chemtype in CARBONS and s.atomname != 'CA'

all_cs = tab.loc[tab.apply(findc, axis=1)]
m, b, r, p, err = linregress(all_cs[EXPT_HEADER], all_cs[CALC_HEADER])
y = lambda x: m*x+b
cmin = min(all_cs[EXPT_HEADER].min(), all_cs[CALC_HEADER].min())
cmax = max(all_cs[EXPT_HEADER].max(), all_cs[CALC_HEADER].max())
mse = ((all_cs[EXPT_HEADER] - all_cs[CALC_HEADER])**2).mean()
ax.text(0.34, 0.08, '$R = %.3f$\n$RMSE = %.2f\ ppm$' % (r, np.sqrt(mse)),
        fontsize=30, transform=ax.transAxes)

ax.plot([cmin, cmax], [y(cmin), y(cmax)], 'k--', linewidth=2)
ax.plot([0, 200], [0, 200], 'k-', linewidth=2)
ax.legend(pls, names, numpoints=1, loc='upper left', prop=legendfont)
ax.set_title('$^{13}$C Shifts', fontdict=dict(fontsize=30))

# Add the nitrogens
ax = fig.add_subplot(224)
ax.tick_params(labelsize=23)
ax.set_title(r'$^{15}$N Shifts', fontdict=dict(fontsize=30))
data = tab.loc[tab.atomname == 'N']
ax.set_xlim([90,150])
ax.set_ylim([90,150])
ax.plot(data[EXPT_HEADER], data[CALC_HEADER], marker='o', color='k',
        markersize=8, ls='None')
ax.plot([90, 150], [90, 150], 'k-', linewidth=2)

m, b, r, p, err = linregress(data[EXPT_HEADER], data[CALC_HEADER])
y = lambda x: m*x+b
hmin = data[EXPT_HEADER].min()
hmax = data[EXPT_HEADER].max()
mse = ((data[EXPT_HEADER] - data[CALC_HEADER])**2).mean()
ax.text(0.34, 0.08, '$R = %.3f$\n$RMSE = %.2f\ ppm$' % (r, np.sqrt(mse)),
        fontsize=30, transform=ax.transAxes)
ax.plot([hmin, hmax], [y(hmin), y(hmax)], 'k--', linewidth=2)

# Plot it
fig.tight_layout()
#plt.show()
fig.savefig(OUTPUT_FIGURE_NAME)
