#!/bin/bash
#rm temp.out #DEBUG
# Solvate.sh: Solvate each frame with target number of waters
# Daniel R. Roe 2015
CPPTRAJ=`which cpptraj`
if [ -z "$CPPTRAJ" ] ; then
  echo "Error: cpptraj not found." > /dev/stderr
  exit 1
fi
TLEAP=`which tleap`
if [ -z "$TLEAP" ] ; then
  echo "Error: tleap not found." > /dev/stderr
  exit 1
fi
if [ -z "$AMBERHOME" ] ; then
  echo "Error: AMBERHOME is not set." > /dev/stderr
  exit 1
fi
# Determine if new leaprc files present
AMBER_VERSION='old'
if [ -f "$AMBERHOME/dat/leap/cmd/leaprc.water.tip3p" ] ; then
  echo "New LEAPRC files present (Amber 16)."
  AMBER_VERSION='new'
fi

target=""               # Target number of waters
buffer="10.0"           # Initial guess for buffer
bufx="10.0"
bufy="10.0"
PDB=""                  # Input PDB name
TOP="solvated.parm7"    # Output topology name
CRD="solvated.rst7"     # Output coordinates name
LEAPIN=""               # Base leap input script
IONSIN=""               # Additional script for adding ions etc
TEMPLEAP="temp.leap.in" # Temporary leap input script
LEAPOUT="temp.leap.out" # Temporary leap output script
tol=2                   # Tolerance (# of waters off from target allowed)
MODE=0                  # 0 - SolvateOct, 1 - SolvateBox, 2 - SolvateBoxXYZ, 3 - SolvateBoxZ
LOADPDB="yes"           # If yes, use 'loadpdb PDB'; otherwise LEAPIN should set up unit MOLNAME
LOADCMD='loadpdb'       # Command to load input file if LOADPDB is yes
CLOSENESS="0.9"         # Leap closeness parameter, default 1.0
MOLNAME="m"             # Molecule name. 'm' by default.
SOLUTERES=""            # Number of SOLUTE residues. If blank try to guess from PDB
WATER_MODEL=0           # Water model to use for solvation, default tip3p
SU=""                   # Solvent unit
# ------------------------------------------------------------------------------
SOLVENT_UNIT[0]="TIP3PBOX"
SOLVENT_UNIT[1]="SPCBOX"
SOLVENT_UNIT[2]="OPCBOX"
SOLVENT_UNIT[3]="TIP4PEWBOX"
SOLVENT_UNIT[4]=""
if [ "$AMBER_VERSION" = 'old' ] ; then
  SOLVENT_CMD[0]=""
  SOLVENT_CMD[1]="WAT=SPC\nloadAmberParams frcmod.spce\n"
  SOLVENT_CMD[2]="WAT=OPC\nloadAmberParams frcmod.opc\n"
  SOLVENT_CMD[3]="WAT=TP4\nloadAmberParams frcmod.tip4pew\n"
else
  SOLVENT_CMD[0]="source leaprc.water.tip3p\n"
  SOLVENT_CMD[1]="source leaprc.water.spce\n"
  SOLVENT_CMD[2]="source leaprc.water.opc\n"
  SOLVENT_CMD[3]="source leaprc.water.tip4pew\n"
fi
# ------------------------------------------------------------------------------
Help() {
  echo "Usage: Solvate.sh <input_file>"
  echo "Input File Options: (default)"
  echo "    target <#>         ) Target # of waters to add."
  echo "    buffer <buf>       ) Initial buffer size (10.0)."
  echo "    bufx   <buf>       ) Initial buffer X size (mode 2|3 only, 10.0)."
  echo "    bufy   <buf>       ) Initial buffer Y size (mode 2|3 only, 10.0)."
  echo "    pdb <file>         ) Solute PDB file name." 
  echo "    top <name>         ) Output topology (solvated.parm7)."
  echo "    crd <name>         ) Output coordinates (solvated.rst7)." 
  echo "    leapin <file>      ) Leap input script for loading parameters etc."
  echo "    ionsin <file>      ) Optional Leap input for loading ions etc (run after solvating)."
  echo "    templeap <name>    ) Name of temporary leap input script (temp.leap.in)." 
  echo "    tol <#>            ) Number of waters > target allowed, will be removed (2)."
  echo "    mode <#>           ) Solvate mode:"
  echo "                         (0)- SolvateOct"
  echo "                          1 - SolvateBox"
  echo "                          2 - SolvateBoxXYZ (bufx and bufy are scaled)"
  echo "                          3 - SolvateBoxZ (bufx and bufy are fixed)"
  echo "    loadpdb {yes|no}   ) If (yes), use 'loadpdb PDB'; otherwise <leapin> should set up unit <molname>."
  echo "    loadcmd <cmd>      ) Command to load solute file; default 'loadpdb'."
  echo "    soluteres <#>      ) Number of solute residues. If blank try to guess from PDB."
  echo "    molname <name>     ) Solute molecule unit name ('m')."
  echo "    solventunit <name> ) Solvent unit (TIP3PBOX)."
  printf "      Recognized solvent units:"
  i=0
  while [[ ! -z ${SOLVENT_UNIT[$i]} ]] ; do
    printf " ${SOLVENT_UNIT[$i]}"
    ((i++))
  done
  echo ""
}
# ------------------------------------------------------------------------------
# ReadInput <file>
ReadInput() {
  if [[ -z $1 || ! -f "$1" ]] ; then
    echo "Error: Solvate Input File '$1' not found." > /dev/stderr
    exit 1
  fi
  # Read options from file
  while read OPTLINE ; do
    #echo "OptionLine: {$OPTLINE}"
    OPT=`echo $OPTLINE | awk '{ print $1; }'`
    # Skip blanks and comments
    if [[ -z $OPT || ${OPT:0:1} = "#" ]] ; then
      continue
    fi
    VAR=`echo $OPTLINE | awk 'BEGIN{ token=1; }{
      for (col=2; col <= NF; col++) {
        if (token > 1) printf(" ");
        if ($col != "") { printf("%s", $col); token++; }
      }
    }'`
    #echo "    Option: $OPT  Variable: $VAR"
    eval VAR=$VAR
    case "$OPT" in
      target      ) target=$VAR ;;
      buffer      ) buffer=$VAR ;;
      bufx        ) bufx=$VAR ;;
      bufy        ) bufy=$VAR ;;
      pdb         ) PDB=$VAR ;;
      top         ) TOP=$VAR ;;
      crd         ) CRD=$VAR ;;
      leapin      ) LEAPIN=$VAR ;;
      ionsin      ) IONSIN=$VAR ;;
      templeap    ) TEMPLEAP=$VAR ;;
      tol         ) tol=$VAR ;;
      mode        ) MODE=$VAR ;;
      loadpdb     ) LOADPDB=$VAR ;;
      loadcmd     ) LOADCMD=$VAR ;;
      soluteres   ) SOLUTERES=$VAR ;;
      molname     ) MOLNAME=$VAR ;;
      solventunit ) SU=$VAR ;;
      * ) echo "Unrecognized options: $OPT" > /dev/stderr ; exit 1 ;;
    esac
  done < $1
  # Check options
  if [[ -z $target ]] ; then
    echo "Error: specify target # of waters 'target <n>'" > /dev/stderr
    exit 1
  fi
  if [[ ! -f "$LEAPIN" ]] ; then
    echo "Error: specify leap script for loading parameters 'leapin <file>'." > /dev/stderr
    exit 1
  fi
  if [[ "$MODE" != "0" && "$MODE" != "1" && "$MODE" != "2" && "$MODE" != "3" ]] ; then
    echo "Error: Specify box mode 'mode <#>': " > /dev/stderr
    echo "Error:   0=SolvateOct, 1=SolvateBox, 2=SolvateBoxXYZ, 3=SolvateBoxZ" > /dev/stderr
    exit 1
  fi
  if [[ $LOADPDB != "yes" && $LOADPDB != "no" ]] ; then
    echo "Error: loadpdb must be 'yes' or 'no'." > /dev/stderr
    exit 1
  fi
  if [[ $LOADPDB = "yes" && -z $PDB ]] ; then
    echo "Error: specify initial non-solvated PDB 'pdb <file>'." > /dev/stderr
    exit 1
  fi
  if [[ -z $SOLUTERES ]] ; then
    if [[ -z $PDB ]] ; then
      echo "Error: If 'pdb <file>' not specified 'soluteres <# solute res>' must be specified." > /dev/stderr
      exit 1
    fi
    if [[ ! -e $PDB ]] ; then
      echo "Error: 'soluteres <# solute res>' not specified and $PDB not found." > /dev/stderr
      exit 1
    fi
  fi 
  # Solvent unit
  if [[ ! -z $SU ]] ; then
    # Try to match up the solvent unit.
    WATER_MODEL=-1
    WM=0
    while [[ ! -z ${SOLVENT_UNIT[$WM]} ]] ; do
      if [[ ${SOLVENT_UNIT[$WM]} = $SU ]] ; then
        echo "Using solvent unit $SU."
        WATER_MODEL=$WM
        break
      fi
      ((WM++))
    done
    if [[ $WATER_MODEL -eq -1 ]] ; then
      echo "Unknown solvent unit: $SU" > /dev/stderr
      exit 1
    fi
  fi
}
# ------------------------------------------------------------------------------
Solvate() {
  if [[ $LOADPDB = "yes" ]] ; then
    echo "$MOLNAME = $LOADCMD $PDB" >> $TEMPLEAP
  fi
  if [[ ! -z ${SOLVENT_CMD[$WATER_MODEL]} ]] ; then
    printf "${SOLVENT_CMD[$WATER_MODEL]}" >> $TEMPLEAP
  fi
  SU=${SOLVENT_UNIT[$WATER_MODEL]}
  case $MODE in 
    0 ) echo "solvateoct $MOLNAME $SU $buffer $CLOSENESS" >> $TEMPLEAP ;;
    1 ) echo "solvatebox $MOLNAME $SU $buffer $CLOSENESS" >> $TEMPLEAP ;;
    2  | 3 )
        echo "solvatebox $MOLNAME $SU {$bufx $bufy $buffer} $CLOSENESS" >> $TEMPLEAP ;;
    * ) echo "Error: Unrecognized mode ($MODE)" > /dev/stderr
        exit 1
    ;;
  esac
  if [[ ! -z $IONSIN && -e $IONSIN ]] ; then
    cat $IONSIN >> $TEMPLEAP
  fi
}
# ------------------------------------------------------------------------------
SaveAndQuit() {
  echo "saveamberparm $MOLNAME $TOP $CRD" >> $TEMPLEAP
  echo "quit" >> $TEMPLEAP
}
# ------------------------------------------------------------------------------
RunTleap() {
  $TLEAP -f $TEMPLEAP > $LEAPOUT
  if [[ $? -ne 0 || ! -f "$TOP" || ! -s "$TOP" ]] ; then
    echo ""
    echo "Error: leap failed. Check '$LEAPOUT' for more details." > /dev/stderr
    exit 1
  fi
}
# ------------------------------------------------------------------------------
CheckBuffer() {
  if [[ ! -z $1 ]] ; then
    if [[ `echo "$1 < 0.00001" | bc` -eq 1 ]] ; then
      echo "Error: New buffer is negative or too small ($1)." > /dev/stderr
      exit 1
    fi
  else
    echo "Error: Buffer is blank." > /dev/stderr
    exit 1
  fi
}
# ------------------------------------------------------------------------------
if [[ $1 = "-h" || $1 = "--help" ]] ; then
  Help
  exit 0
fi
ReadInput $1
loop=1        # Keep guessing while loop==1
ntries=0      # To keep track of how many tries
lastdiff=0   
# Get number of solute residues in PDB
if [[ -z $SOLUTERES ]] ; then
  SOLUTERES=`$CPPTRAJ -p $PDB -mr '*' | awk '{print NF - 1;}'`
  if [[ -z $SOLUTERES ]] ; then
    echo "Error: Could not determine # solute residues." > /dev/stderr
    exit 1
  fi
fi
echo "Solute has $SOLUTERES residues"
 
while [[ $loop -eq 1 ]] ; do
  if [[ $MODE -ge 2 ]] ; then
    printf "%i) Buffer: %f %f %f" $ntries $bufx $bufy $buffer
  else
    printf "%i) Buffer: %f" $ntries $buffer
  fi
  # Set up buffer
  cp $LEAPIN $TEMPLEAP
  Solvate
  SaveAndQuit
  # Run tleap with given buffer
  RunTleap
  # Check that a coordinate file was actually written
  if [[ ! -f "$CRD" || ! -f "$TOP" ]] ; then
    echo "LEAP error: coordinate/topology files not written." > /dev/stderr
    exit 1
  fi
  # How many waters were added?
  solvent_res_added=`$CPPTRAJ -p $TOP -mr ':WAT' | awk '{print NF - 1}'`
  if [[ -z $solvent_res_added ]] ; then
    echo "Error: Could not get number of solvent residues." > /dev/stderr
    exit 1
  fi
  # How far off is it?
  ((diff = target - solvent_res_added))
  printf " #Res: %i  Diff: %i " $solvent_res_added $diff
  # If this is the first time through choose an appropriate change val
  if [[ $ntries -eq 0 ]] ; then
    #change=$(echo "$diff / 100" | bc -l)
    change="0.001"
  fi
  # See if we have tol more waters than the target.
  negtol=0
  if [[ $tol -gt 0 ]] ; then
    ((negtol = -$tol))
  fi
  if [[ $diff -lt 0 && $diff -ge $negtol ]] ; then
    # Close enough, just tell leap to remove the last water residue(s)
    ((LASTRES = $SOLUTERES + $solvent_res_added))
    ((FIRSTRES = $LASTRES + $diff + 1))
    echo "Only $diff off, removing solvent residue(s) ($FIRSTRES-$LASTRES)"
    cp $LEAPIN $TEMPLEAP
    Solvate
    for ((RRES = $FIRSTRES; RRES <= $LASTRES; RRES++)) ; do
      echo "remove $MOLNAME $MOLNAME.$RRES" >> $TEMPLEAP
    done
    SaveAndQuit
    RunTleap
    echo "$solvent_res_added      $buffer      $diff" >> buffers.dat
    exit 0
  elif [[ $diff -eq 0 ]] ; then
    # If diff==0 then we're done
    loop=0
    echo -n "Found"
    echo "$solvent_res_added      $buffer" >> buffers.dat
  else
    #echo $diff >> temp.out #DEBUG
    # If we passed the zero point reduce change value.
    CHANGE_DIR="0"
    if [[ $lastdiff -gt 0 && $diff -lt 0 ]] ; then
      change=`echo "$change / 2" | bc -l`
      CHANGE_DIR="-"
    elif [[ $lastdiff -lt 0 && $diff -gt 0 ]] ; then
      change=`echo "$change / 2" | bc -l`
      CHANGE_DIR="-"
    elif [[ $ntries -gt 0 ]] ; then
       # If we took a step and the size of the step we just took is smaller
       # than the remaining distance to target, increase the step size.
       ((STEPSIZE = $diff - $lastdiff))
       if [[ $STEPSIZE -lt 0 ]] ; then
         ((STEPSIZE = -$STEPSIZE))
       fi
       ABSDIFF=$diff
       if [[ $ABSDIFF -lt 0 ]] ; then
         ((ABSDIFF = -$ABSDIFF))
       fi
       #printf " STEP=%i ABS=%i " $STEPSIZE $ABSDIFF
       if [[ $STEPSIZE -gt 0 && $STEPSIZE -lt $ABSDIFF ]] ; then
        change=`echo "$change * 2.4" | bc -l`
        CHANGE_DIR="+"
       fi
    fi
    lastdiff=$diff
    printf "Change: %G ($CHANGE_DIR)  " $change
    # Choose a new buffer value
    newbuffer=`echo "$buffer + ($change * $diff)" | bc -l` # TODO should be minus
    CheckBuffer $newbuffer
    lastbuffer=$buffer
    buffer=$newbuffer
    if [[ $MODE -eq 2 ]] ; then
      newbufx=`echo "$bufx + ($change * $diff)" | bc -l`
      CheckBuffer $newbufx
      newbufy=`echo "$bufy + ($change * $diff)" | bc -l`
      CheckBuffer $newbufy
      bufx=$newbufx
      bufy=$newbufy
    fi
  fi
  # Safety valve
  if [[ $ntries -gt 100 ]] ; then
    echo "Taking too long - moving on..."
    rm $CRD
    loop=0
    exit 1
  fi
  ((ntries++))
  echo ""
  #break # DEBUG
done
exit 0
