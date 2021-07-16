
include Makefile.in

OBJ_DIR = obj

MYEXEC = readGeom.exe

MYFUN = modGzip.f90     \
        readGeom.f90 

MYOBJ = $(patsubst %.f90,$(OBJ_DIR)/%.o,$(MYFUN))


#---- Rules ----#

.PHONY: $(MYEXEC)
$(MYEXEC): $(MYOBJ)
	$(FORTRAN) $^ $(LFLAGS) -o $@

$(MYOBJ): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: %.f90
	$(FORTRAN) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(MYEXEC) $(OBJ_DIR) *.o *.mod
