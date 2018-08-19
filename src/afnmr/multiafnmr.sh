#!/bin/sh

# All shifts-related stuff should be in the directory containing this script
BINDIR=`dirname $0`

#  MacOS doesn't have a "seq" command:
if [ -x /usr/bin/seq ]; then
	seq="/usr/bin/seq"
elif [ -x /usr/bin/jot ]; then
	seq="/usr/bin/jot"
else
	echo "Cannot find either /usr/bin/seq or /usr/bin/jot; exiting"
	exit 1
fi

usage() {

   test $# -gt 0 && echo $1 && echo ""
   echo "`basename $0` <PDB_Basename> [afnmr flags]"
   echo "   This script takes an input PDB file (PDB_Basename.pdb) and "
   echo "   splits it into a separate PDB file for each model (structure)"
   echo "   present in the original. If additional arguments exist, afnmr"
   echo "   will be run with those arguments."
   echo ""
   echo "[AFNMR flags]"
   $BINDIR/afnmr --help | awk '(NR > 6)'
   exit 1
}

error() {
   echo $1; exit 1
}

test "$1" = "-h" && usage
test "$1" = "--help" && usage
test "$1" = "-H" && usage

# Make sure we got an argument
test $# -gt 0 || usage "No arguments found"

BASENAME=$1
shift;

# Make sure our PDB exists
test -f $BASENAME.pdb || usage "$BASENAME.pdb not found"

test -x `which cpptraj` || error "Could not find cpptraj. Install AmberTools."

# Determine how many models were present and split by cpptraj
nframes=`cpptraj -p $BASENAME.pdb -y $BASENAME.pdb -tl | awk '{print $2}'`

cpptraj -p $BASENAME.pdb << EOF
trajin $BASENAME.pdb
trajout $BASENAME.pdb multi
EOF

echo "Found $nframes structures in $BASENAME.pdb."

for i in `$seq $nframes`; do
   /bin/mv $BASENAME.pdb.$i ${BASENAME}${i}.pdb
   echo "  ${BASENAME}${i}.pdb"
done

# If we have no more arguments, bail here
test $# -eq 0 && exit 0

# Run afnmr with remaining arguments in parallel
test -z "$NCPUS" && NCPUS=1

ARGS="$@ -workdir"
python -c "print(' '.join(['${BASENAME}%d' % i for i in range(1,$nframes+1)]))" | \
   xargs -n 1 -P $NCPUS ${BINDIR}/afnmr $ARGS
