#!/bin/bash

# Initial variables... adjust as needed
BMRB=5371
SYS=1LC6
NSTRUCT=40
NRES=24

fake_rdb() {
# This fakes a RDB file with filler numbers if the output file does not exist
   echo "Faking a RDB for \"$1\"" > /dev/stdout
   sys=`python -c "print '$1'[:-3]"`
   resnum=`python -c "print '%d' % (int('$1'[-3:]))"`
   awk "\$5==\"$resnum\" {printf(\"$sys\t$resnum\t%s\t%s\t%f\n\",\$3,\$4,-99999)}" $1.pqr >> $1.rdb
}

# Format nres for zeroes
fNRES=`printf "%03d" $NRES`

echo "Sorting and separating experimental shifts..."
sed -e "s/'/xX/g" bmr${BMRB}.rdb | sorttbl -n res atomname | sed -e "s/xX/'/g" | \
	row atomname mat /H/ > bmr${BMRB}_HNMR.rdb
sed -e "s/'/xX/g" bmr${BMRB}.rdb | sorttbl -n res atomname | sed -e "s/xX/'/g" | \
	row atomname mat /C/ > bmr${BMRB}_CNMR.rdb
sed -e "s/'/xX/g" bmr${BMRB}.rdb | sorttbl -n res atomname | sed -e "s/xX/'/g" | \
	row atomname mat /N/ > bmr${BMRB}_NNMR.rdb

echo "Extracting shifts..."

namelist=""
for i in `eval echo "{1..$NSTRUCT}"`; do
   for j in `eval echo "{001..$fNRES}"`; do
      namelist="$namelist ${SYS}${i}${j}"
   done
done

for basename in $namelist; do
   echo "sys	res	atomname	resname	demon3-shift" > $basename.rdb
   echo "" >> $basename.rdb
   if [ -f $basename.pqr -a -f $basename.out ]; then
      getshifts-demon3 $basename | sed -e "s/ //g" >> $basename.rdb
   else
      fake_rdb $basename
   fi
done

echo "Combining all residues and separating nuclei..."

for i in `eval echo "{1..$NSTRUCT}"`; do
   cat ${SYS}${i}001.rdb > tmp
   for j in `eval echo "{002..$fNRES}"`; do
      headchg -del < ${SYS}${i}${j}.rdb >> tmp
   done
   sed -e "s/'/xX/g" tmp | sorttbl -n res atomname | row atomname mat /H/ | \
         column res atomname resname demon3-shift | \
         sed -e "s/xX/'/g" > ${SYS}${i}_HNMR.rdb
   sed -e "s/'/xX/g" tmp | sorttbl -n res atomname | row atomname mat /C/ | \
         column res atomname resname demon3-shift | \
         sed -e "s/xX/'/g" > ${SYS}${i}_CNMR.rdb
   sed -e "s/'/xX/g" tmp | sorttbl -n res atomname | row atomname mat /N/ | \
         column res atomname resname demon3-shift | \
         sed -e "s/xX/'/g" > ${SYS}${i}_NNMR.rdb

   /bin/rm -f tmp
done

echo "Joining all H-NMR frames..."

cmd="jointbl < ${SYS}2_HNMR.rdb res atomname resname ${SYS}1_HNMR.rdb"
for i in `eval echo "{3..$NSTRUCT}"`; do
   cmd="$cmd | jointbl res atomname resname ${SYS}${i}_HNMR.rdb"
done
cmd="$cmd | column -v res atomname resname | headchg -del > combined_HNMR.rdb"
eval $cmd

echo "Joining all C-NMR frames..."

cmd="jointbl < ${SYS}2_CNMR.rdb res atomname resname ${SYS}1_CNMR.rdb"
for i in `eval echo "{3..$NSTRUCT}"`; do
   cmd="$cmd | jointbl res atomname resname ${SYS}${i}_CNMR.rdb"
done
cmd="$cmd | column -v res atomname resname | headchg -del > combined_CNMR.rdb"
eval $cmd

echo "Joining all N-NMR frames..."

cmd="jointbl < ${SYS}2_NNMR.rdb res atomname resname ${SYS}1_NNMR.rdb"
for i in `eval echo "{3..$NSTRUCT}"`; do
   cmd="$cmd | jointbl res atomname resname ${SYS}${i}_NNMR.rdb"
done
cmd="$cmd | column -v res atomname resname | headchg -del > combined_NNMR.rdb"
eval $cmd

echo "Computing shift average/variance"

python << EOF
import numpy as np
from rdb import RDB

def _average(vector):
   """ Skip -99999 fillers in average """
   return np.asarray([i for i in vector if i != -99999]).mean()
def _stdev(vector):
   """ Skip -99999 fillers in average """
   return np.asarray([i for i in vector if i != -99999]).std()

with open('${SYS}1_HNMR.rdb', 'r') as f:
   f.readline(); f.readline()
   hatname = ['%s\\t%s' % tuple(line.split()[0:2]) for line in f]
with open('${SYS}1_CNMR.rdb', 'r') as f:
   f.readline(); f.readline()
   catname = ['%s\\t%s' % tuple(line.split()[0:2]) for line in f]
with open('${SYS}1_NNMR.rdb', 'r') as f:
   f.readline(); f.readline()
   natname = ['%s\\t%s' % tuple(line.split()[0:2]) for line in f]

hnmrdata = np.loadtxt('combined_HNMR.rdb')
cnmrdata = np.loadtxt('combined_CNMR.rdb')
nnmrdata = np.loadtxt('combined_NNMR.rdb')

hdata = RDB()
cdata = RDB()
ndata = RDB()

hdata['sys'] = ['${SYS}' for i in hatname]
hdata['res'] = [int(n.split()[0]) for n in hatname]
hdata['atomname'] = [n.split()[1] for n in hatname]
hdata['qmshift_avg'] = [_average(vec) for vec in hnmrdata]
hdata['qmshift_std'] = [_stdev(vec) for vec in hnmrdata]
hdata.check()
hdata.write_to_file('_combined_HNMR.rdb')

cdata['sys'] = ['${SYS}' for i in catname]
cdata['res'] = [int(n.split()[0]) for n in catname]
cdata['atomname'] = [n.split()[1] for n in catname]
cdata['qmshift_avg'] = [_average(vec) for vec in cnmrdata]
cdata['qmshift_std'] = [_stdev(vec) for vec in cnmrdata]
cdata.write_to_file('_combined_CNMR.rdb')

ndata['sys'] = ['${SYS}' for i in natname]
ndata['res'] = [int(n.split()[0]) for n in natname]
ndata['atomname'] = [n.split()[1] for n in natname]
ndata['qmshift_avg'] = [_average(vec) for vec in nnmrdata]
ndata['qmshift_std'] = [_stdev(vec) for vec in nnmrdata]
ndata.write_to_file('_combined_NNMR.rdb')

EOF

echo "Combining shifts with experiment..."

grep -v "HO[12]'" _combined_HNMR.rdb | jointbl res atomname bmr${BMRB}_HNMR.rdb | \
   column res resname atomname qmshift_avg qmshift_std obs > ${SYS}_HNMR.rdb
jointbl < bmr${BMRB}_CNMR.rdb res atomname _combined_CNMR.rdb | \
   column res resname atomname qmshift_avg qmshift_std obs > ${SYS}_CNMR.rdb
jointbl < bmr${BMRB}_NNMR.rdb res atomname _combined_NNMR.rdb | \
   column res resname atomname qmshift_avg qmshift_std obs > ${SYS}_NNMR.rdb

echo "Done. To clean:"
echo "	/bin/rm -f _combined_?NMR.rdb combined_?NMR.rdb ${SYS}[0-9]*.rdb ${SYS}[1-9][1-9]*_?NMR.rdb"
