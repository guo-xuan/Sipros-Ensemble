#!/bin/bash 

# Get path of run script
if [ -L $0 ] ; then
    exePath=$(dirname $(readlink -f $0)) ;
else
    exePath=$(dirname $0) ;
fi;

# Generate PSM table

# PSM Filtering

# Protein Assembly
