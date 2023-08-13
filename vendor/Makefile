F77 = f77
FFLAGS = -O3
AR = ar

OUT_DIR ?= .

FFILES = lufact.f lusolv.f ludfs.f lucomp.f \
		lucopy.f lsolve.f usolve.f ifill.f \
		rfill.f icopy.f rcopy.f requiv.f \
		second.f flops.f \
		ltsolve.f utsolve.f maxmatch.f ordstat.f \
		cntrow.f

OFILES = $(addprefix $(OUT_DIR)/, $(addsuffix .o,$(basename $(FFILES))))

LIB = $(OUT_DIR)/libgp.a

all: $(LIB)

$(LIB): $(OFILES)
	$(AR) cr $@ $(OFILES)

$(OUT_DIR)/%.o: %.f
	$(F77) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OFILES)

purge: clean
	rm -f $(LIB)
