#!/bin/bash 

#Get path of run script
if [ -L $0 ] ; then
    exePath=$(dirname $(readlink -f $0)) ;
else
    exePath=$(dirname $0) ;
fi;

dataOutPath=""
workingPath=""
ms2Path=""
configGroupPath=""
configPath=""

while [[ $# > 0 ]]
do
key="$1"
case $key in
    -h|--help)              # help output
    echo -e "Usage:\n"
    echo -e "   runSipros.sh [OPTION]...<PARAM>...\n\n"
    echo -e "<PARAMS>\n"
    echo -e "   -w\t directory name with multiple MS2 files (required if -f not set).\n"
    echo -e "   -f\t single MS2 filename (required if -w not set).\n"
    echo -e "   -g\t directory name with multiple configuration files (required -f -c not set).\n"
    echo -e "   -c\t configuration filename (required if -g not set).\n"
    echo -e "   -o\t output directory path (required).\n"
    echo -e "<OPTIONS>\n"
    echo -e "   -h\t help.\n"
    exit 1
    ;;
    -o)						# output directory
    dataOutPath="$2"
    shift # past argument
    ;;
    -w)						# working directory
    workingPath="$2"
    shift # past argument
    ;;
    -f)						# single ms2 file
    ms2Path="$2"
    shift # past argument
    ;;
    -g)						# configuration directory
    configGroupPath="$2"
    shift # past argument
    ;;
    -c)						# single configuration file
    configPath="$2"
    shift # past argument
    ;;
    *)
    echo "ERROR: Unidentified user variable $key"
    exit 1        				# unknown option
    ;;
esac
shift # past argument or value
done

if [ -z "$ms2Path" ] && [ -z "$workingPath" ] ; then
   echo "Input incomplete. Not all required parameters specified. Run with -h for help. Exiting..."
   exit 1
fi

if [ -z "$configGroupPath" ] && [ -z "$configPath" ] ; then
   echo "Input incomplete. Not all required parameters specified. Run with -h for help. Exiting..."
   exit 1
fi

if [ -z "$dataOutPath" ] ; then
   echo "Input incomplete. Not all required parameters specified. Run with -h for help. Exiting..."
   exit 1
fi

#Check required binaries are in the bin directory
SIPROSOPENMP="${exePath}/bin/Sipros_OpenMP"
if [ -f $SIPROSOPENMP ] ; then
   echo "$SIPROSOPENMP exists."
else
   echo "$SIPROSOPENMP does not exist in script directory."
fi

SIPROSMPI="${exePath}/bin/Sipros_MPI"
if [ -f $SIPROSMPI ] ; then
   echo "$SIPROSMPI exists."
else
   echo "$SIPROSMPI does not exist in script directory."
fi

#check if output directory exists, if not create it
if [ -d $dataOutPath ] ; then
   echo "Output directory exists."
else
   echo "Cresting output directory: $dataOutPath"
   `mkdir $dataOutPath`
fi

logFile="${dataOutPath}/$(date).log"

echo Starting Time is $(date)

#database searching
if [ -z "$ms2Path" ] ; then
   if [ -z "$configGroupPath" ] ; then
      ${exePath}/bin/Sipros_OpenMP -w ${workingPath} -c ${configPath} -o ${dataOutPath} -s > ${logFile}
   else
      ${exePath}/bin/Sipros_OpenMP -w ${workingPath} -g ${configGroupPath} -o ${dataOutPath} -s > ${logFile}
   fi
else
   if [ -z "$configGroupPath" ] ; then
      ${exePath}/bin/Sipros_OpenMP -f ${ms2Path} -c ${configPath} -o ${dataOutPath} -s > ${logFile}
   else
      ${exePath}/bin/Sipros_OpenMP -f ${ms2Path} -g ${configGroupPath} -o ${dataOutPath} -s > ${logFile}
   fi
fi

echo Ending Time is $(date)