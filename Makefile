SRCS = src/ParametricMatrixTools.mpl
MLAFILE = ParametricMatrixTools.mla
ASSERTLEVEL = 2

all: mla

mint: 
	@$(foreach f,$(SRCS), mint -q $(MACROS) $(f);)
	@(echo "Done!";)

test: mla
	@maple -q -A $(ASSERTLEVEL) -B -b $(MLAFILE) Tests/tests.tst

mla: $(SRCS)
	@rm -f $(MLAFILE)
	@maple -q src/make_mla.mpl
	@(echo "Done!";)

clean:
	rm -f $(MLAFILE)