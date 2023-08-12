F77 = f77
FFLAGS = -O1

AR = ar

LUFILES = lufact.o lusolv.o ludfs.o lucomp.o \
		lucopy.o lsolve.o usolve.o ifill.o \
		rfill.o icopy.o rcopy.o requiv.o \
		second.o flops.o \
                ltsolve.o utsolve.o maxmatch.o ordstat.o \
                cntrow.o

LIB = libgp.a

all: $(LIB)

$(LIB): $(LUFILES)
	$(AR) cr $@ $(LUFILES)

%.o: %.f
	$(F77) $(FFLAGS) -c $<

install: $(LIB)
	cp $< $(OUT_DIR)

clean:
	rm -f $(LUFILES)

purge: clean
	rm -f $(LIB)
