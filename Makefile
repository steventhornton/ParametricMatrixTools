SRCS = ParametricMatrixTools.mpl
MLAFILE = ParametricMatrixTools.mla
MACROS = -D RC=RegularChains\
		 -D RC_CST=RegularChains:-ConstructibleSetTools\
		 -D RC_PST=RegularChains:-ParametricSystemTools\
		 -D RC_CT=RegularChains:-ChainTools\
		 -D RC_SAST=RegularChains:-SemiAlgebraicSetTools\
		 -D RC_MT=RegularChains:-MatrixTools\
		 -D LA=LinearAlgebra\
		 -D PT=PolynomialTools\
		 -D LT=ListTools
ASSERTLEVEL = 2

all: mla

mint: 
	@$(foreach f,$(SRCS), mint -q $(MACROS) $(f);)
	@(echo "Done!";)

test: mla
	@maple -q -A $(ASSERTLEVEL) -B -b $(MLAFILE) Tests/tests.tst

mla: $(SRCS)
	@rm -f $(MLAFILE)
	@maple -q $(MACROS) make_mla.mpl
	@(echo "Done!";)

clean:
	rm -f $(MLAFILE)