
include Makefile.in

MYEXEC = myzlib.exe

MYFUN = myzlib_f.f90     \
        myzlib_c.cpp 

MYOBJ = myzlib_f.o     \
        myzlib_c.o

#---- Rules ----#

.PHONY: $(MYEXEC)
$(MYEXEC): $(MYOBJ)
	$(FORTRAN) $^ $(LFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.f90
	$(FORTRAN) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(MYEXE) *.o *.mod
