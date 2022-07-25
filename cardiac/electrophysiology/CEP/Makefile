
.SUFFIXES: .f .o

CUR_DIR:= $(dir $(abspath $(lastword $(MAKEFILE_LIST)) ))
SRC_DIR = $(CUR_DIR)src
OBJ_DIR = $(CUR_DIR)obj
BIN_DIR = $(CUR_DIR)bin

FC = gfortran
FFLAGS  = -O3
FFLAGS += -J $(OBJ_DIR)
LAPACK_LFLAGS = -llapack

debug = 0
ifeq ($(debug),1)
   FFLAGS += -g -fbacktrace -fbounds-check -Wall
endif

INCLUDES += -I./include

MYFUN  = MATFUN.f    \
         UTIL.f      \
         APMOD.f     \
         FNMOD.f     \
         TTPMOD.f    \
         BOMOD.f     \
         MOD.f       \
         MAIN.f

MYEXE  = $(BIN_DIR)/cep.exe

MYSRC = $(patsubst %,   $(SRC_DIR)/%,   $(MYFUN))
MYOBJ = $(patsubst %.f, $(OBJ_DIR)/%.o, $(MYFUN))

OPTIONS = $(FFLAGS)  $(LAPACK_LFLAGS)

#---Rules---#
.PHONY: $(MYEXE)
$(MYEXE): $(MYOBJ)  $(BIN_DIR)
	$(FC) $(MYOBJ) $(OPTIONS) -o $@

$(MYOBJ): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR)

cleanall:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

