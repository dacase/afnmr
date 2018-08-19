#!/usr/bin/env python

from __future__ import print_function, with_statement, division

from argparse import ArgumentParser
try:
   from collections import OrderedDict
except ImportError:
   from ordereddict import OrderedDict
import os
import re
import sys as _sys

parser = ArgumentParser()

parser.add_argument('basename', help='''Name of the system to analyze.
            <basename>.emp, <basename>.pdb, and <basename>###.pqr must exist,
            where ### is the residue number from 1 to the number of
            residues.''')
parser.add_argument('--verbose', help='''Print out all skipped ring current
            contributions (since they were detected in the QM region''',
            action='store_true', default=False)
parser.add_argument('-f', '--scaling-factor', default=0.0, type=float,
            help='''The fraction of the ring current contribution for rings
            \033[1minside\033[m the quantum region. Default is %(default)f''',
            metavar='FLOAT')
parser.add_argument('-c', '--coil', default=False, action='store_true',
            help='''Print the random coil results to the resulting RDB file
            (which is printed to stdout)''')
parser.add_argument('-s', '--split-qm-corrections', default=False, dest='split',
            action='store_true', help='''Print the corrections from the
            classical region and the quantum region to separate columns.''')

opt = parser.parse_args()

VERBOSE = opt.verbose

# If we asked to split QM corrections and scaling factor is 0, set the scaling
# factor to 1 to print out the full QM correction so it can be later scaled via
# post-processing
if opt.split and opt.scaling_factor == 0:
   opt.scaling_factor = 1.0

# Add AMBERHOME/bin to the PYTHONPATH
_sys.path.append(os.path.join(os.getenv('AMBERHOME'), 'bin'))

from chemistry import system

fullsys = system.ChemicalSystem.load_from_pdb('%s.pdb' % opt.basename)

nres = 0
for res in fullsys:
   if res.name in ('HOH', 'WAT'): continue
   nres += 1

# Now store all of the systems for all of the residues
syslist = [system.ChemicalSystem.load_from_pqr('%s%03d.pqr' % (opt.basename, i))
               for i in range(1, nres+1)]

# Now get the detailed shifts
ret = os.system('shifts -details -reslib "::H*" %s > /dev/null 2>&1' % 
                opt.basename)

if ret != 0:
   _sys.exit('Running shifts failed!')

# Now get the detailed ring contributions
ringcont = OrderedDict()
qmringcont = OrderedDict()
coil = dict() # Ordering taken from ringcont OrderedDict
with open('%s.emp' % opt.basename, 'r') as femp:
   contre = re.compile(r'Detailed contributions for atom '
         r'''([A-Za-z0-9\-'"]+):(\d+):([A-Za-z0-9\-'"]+)''')
   rcre = re.compile(
         r'(\w+) *' # system name
         r'(\w+) *' # residue name
         r'(\d+) *' # residue number
         r'''([\w'"]+) *''' # atom name
         r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?) *' # ring current
         r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?) *' # electrostatic
         r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?) *' # anisotropy
         r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?) *' # constant
         r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?) *' # random coil
         r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?) *' # predicted
   )
   line = femp.readline()
   while line:
      rematch = contre.match(line)
      if not rematch:
         # See if this line has a random coil contribution
         rematch2 = rcre.match(line)
         if rematch2:
            gr = rematch2.groups()
            curkey = '%d %s' % (int(gr[2]), gr[3])
            coil[curkey] = float(gr[8])
         # Skip to the next line and continue
         line = femp.readline()
         continue
      femp.readline() # eat a blank line
      resname, resnum, atomname = rematch.groups()
      resnum = int(resnum)
      sys = syslist[resnum - 1]
      qm_included_residues = [res.number for res in sys if len(res) > 1]
      curkey = '%d %s' % (resnum, atomname)
      ringcont[curkey] = 0.0
      qmringcont[curkey] = 0.0
      line = femp.readline()
      while line[:5] == 'ring:':
         if not line:
            sys.exit('Unexpected EOF when parsing!')
         words = line.split()
         resnum, cont = int(words[2]), float(words[3])
         if resnum not in qm_included_residues:
            ringcont[curkey] += cont
         else:
            if VERBOSE:
               print('Detected ring in quantum section for atom %s. Adding %s' %
                     (curkey, opt.scaling_factor * cont),
                     file=_sys.stderr)
            qmringcont[curkey] += opt.scaling_factor * cont
         line = femp.readline()
      line = femp.readline()

# Now we have the contributions
print('# Ring current contributions in QM region scaled by %s' %
      opt.scaling_factor)
if not opt.coil:
   if not opt.split:
      print('res\tatomname\tcorrection\n', file=_sys.stdout)
   else:
      print('res\tatomname\tclassical_correction\tqm_correction\n',
            file=_sys.stdout)
else:
   if not opt.split:
      print('res\tatomname\tcoil\tcorrection\n', file=_sys.stdout)
   else:
      print('res\tatomname\tcoil\tclassical_correction\tqm_correction\n',
            file=_sys.stdout)
   # Fill in missing coil contributions (these are really just holders)
   for key in ringcont:
      if not key in coil: coil[key] = 0.0
for key in ringcont:
   resnum, atomname = key.split()
   if not opt.coil:
      if not opt.split:
         rc = ringcont[key] + qmringcont[key]
         print('%s\t%s\t%s' % (resnum, atomname, rc), file=_sys.stdout)
      else:
         print('%s\t%s\t%s\t%s' % (resnum, atomname, ringcont[key],
            qmringcont[key]), file=_sys.stdout)
   else:
      if not opt.split:
         rc = ringcont[key] + qmringcont[key]
         print('%s\t%s\t%s\t%s' % (resnum, atomname, coil[key], rc),
               file=_sys.stdout)
      else:
         print('%s\t%s\t%s\t%s\t%s' % (resnum, atomname, coil[key],
            ringcont[key], qmringcont[key]), file=_sys.stdout)
