
include Makefile.in

SRC_DIR = src
OBJ_DIR = obj

FFLAGS  += -J $(OBJ_DIR)

MYEXEC = extrude.exe

MYF90FUN = vtkXMLParser.f90   \
           extrude.f90

MYCFUN = zlibDataReader.c 

MYSRC = $(patsubst %,$(SRC_DIR)/%,$(MYF90FUN))
MYOBJ = $(patsubst %.f90,$(OBJ_DIR)/%.o,$(MYF90FUN))
MYOBJ += $(patsubst %.c,$(OBJ_DIR)/%.o,$(MYCFUN))

#---- Rules ----#

.PHONY: $(MYEXEC)
$(MYEXEC): $(MYOBJ) 
	$(FORTRAN) $^ $(LFLAGS) -o $@

$(MYOBJ): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FORTRAN) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(MYEXEC)

