## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += MHDAM_GNUmakefile

# the base name(s) of the application(s) in this directory
ebase = MHDAM

# the location of Chombo lib dir
CHOMBO_HOME =$(HOME)/MS_FLUKSS-3.2.4-HE/Chombo-3.2.4/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRTimeDependent AMRElliptic AMRTools BoxTools 

# relative paths to source code directories
base_dir = .
src_dirs = ../src ../srcIdeal ../problems ../swlism ../mc ../fftpack 
#src_dirs = ../src ../srcIdeal ../problems ../mc ../main ../fftpack 

# input file for 'run' target
INPUT = swlism.inputs


# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# Note: tell Make.rules to not use -q for 'run' target
_space = # # NOTE: keep the space after the "=" or bad things happen
VERBOSE = $(_space)

# application-specific targets


