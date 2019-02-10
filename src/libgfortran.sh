#!/bin/sh

#  find location of libgfortran:

for x in `gfortran -print-search-dirs | grep libraries | \
          sed -e "s/libraries: =//g" -e "s/:/ /g"`; do
    test -f $x/libgfortran.dylib && break
done

echo "-L$x"
