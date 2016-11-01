# ======================================================================= #
# ======================================================================= #
#
# VanishingConstructibleSet.mpl
#
# AUTHOR .... Steven E. Thornton
#                Under the supervision of 
#                Robert M. Corless & Marc Moreno Maza
# EMAIL ..... sthornt7@uwo.ca
# UPDATED ... Jan. 20/2016
#
# Generate the constructible set such that a polynomial vanishes for all 
# points in the constructible set.
#
# CALLING SEQUENCE
#   VanishingConstructibleSet(p, R, options)
#   VanishingConstructibleSet(p, v, R, options)
#   VanishingConstructibleSet(p, cs, R, options)    # To Do
#   VanishingConstructibleSet(p, v, cs, R, options) # To Do
#
# INPUT:
#   p .... Polynomial
#   v .... Variable
#   R .... Polynomial ring
#
# OPTIONS
#   lazy ......... true (default):
#                      - Computation is done in the sense of Kalkbrener
#                        (see RegularChains)
#                   false:
#                      - Computation is done in the sense of Lazard 
#                        (see RegularChains)
# OUTPUT
#    TO DO:
#        If a constructible set is provided as one of the input arguents then the 
#        output constructible set will be a subset of the input constructible set 
#        such that p vanishes at all points in the output constructible set.
#
#   3 Arguments:
#       A constructible set such that all coefficients of p in v are zero.
#       If any of the coefficients are constant then the empty constructible set
#       is returned.
#   2 Arguments:
#       A constructible set such that p is zero. If p is constant and non-zero 
#        then the empty constructible set is returned.
# ======================================================================= #
# ======================================================================= #
VanishingConstructibleSet := module()

    export ModuleApply;

    local
        init,
        init_no_var,
        init_with_var,
        processOptions,
        checkInput;

    ModuleApply := proc() :: TRDcs;
        return init(args);
    end proc:


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ------------------------------------------------------------------------------
# init
#
# Checks the types of the input and calls the implementation if all input values
# pass checks.
#
# INPUT/OUTPUT
#   Same as VanishingConstructibleSet
# ------------------------------------------------------------------------------
init := proc() :: TRDcs;
    
    local p, v, R, opts;
    
    if nargs < 2 then
        error "Insufficient number of arguments";
    elif nargs > 4 then
        error "To many arguments";
    end if;
    
    if type(args[2], 'TRDring') then
        # VanishingConstructibleSet(p, R, options)
        
        p := args[1];
        R := args[2];
        
        opts := processOptions({args[3..-1]});
        
        return init_no_var(p, R, opts);
        
    elif type(args[2], name) then
        # VanishingConstructibleSet(p, v, R, options)
        if nargs < 3 then
            error "Expected third argument of a polynomial ring";
        end if;
        
        p := args[1];
        v := args[2];
        R := args[3];
        
        opts := processOptions({args[4..-1]});
        
        return init_with_var(p, v, R, opts);
        
    else
        error "Expected second argument to be a polynomial ring or a variable";
    end if;
    
end proc:


# ------------------------------------------------------------------------------
# processOptions
#
# Extracts the options from a set, if an option is missing the default values is
# returned.
#
# INPUT
#   A set of equations corresponding to the options.
#
# OUTPUT
#    A table with indices
#        'lazy'
#    See VanishingConstructibleSet header for specifications.
# ------------------------------------------------------------------------------
processOptions := proc(opts_in::set(equation), $) :: table;

    local opts::table(),
          opt::equation;

    # Default values
    opts['lazy'] := true;

    # Process each option
    for opt in opts_in do
        if lhs(opt) in {indices(opts, 'nolist')} then
            opts[lhs(opt)] := rhs(opt);
        else
            error "'%1' is not a valid option", lhs(opt);
        end if;
    end do;

    return opts;

end proc:


# ------------------------------------------------------------------------------
# checkInput
#
# Checks for errors in the input values
#
# INPUT
#   p ...... Polynomial
#   R ...... Polynomial Ring
#    opts ... A table containing the options (see VanishingConstructibleSet 
#             header)
# ------------------------------------------------------------------------------
checkInput := proc(p::polynom, R::TRDring, opts::table, $)

    # p must be a polynomial in R
    if not RC:-TRDis_poly(p, R) then
        error "invalid polynomial: %1", p;
    end if;
    
    # Check the lazy option
    if not type(opts['lazy'], 'truefalse') then
        error "lazy option must be a boolean values";
    end if;
    
end proc:


# ------------------------------------------------------------------------------
# init_no_var
#
# Compute a triangular decomposition of p such that p vanishes at every
# point in the triangular decomposition.
#
# INPUT
#   p ...... Polynomial
#   R ...... Polynomial Ring
#    opts ... A table containing the options (see VanishingConstructibleSet 
#             header)
#
# OUTPUT
#    Same as VanishingConstructibleSet
# ------------------------------------------------------------------------------
init_no_var := proc(p::polynom, R::TRDring, opts::table, $) :: TRDcs;
    
    local lrc::list(TRDrc), lrs::TRDlrs;
    
    # Check the input for errors
    checkInput(p, R, opts);
    
    if opts['lazy'] then
        lrc := RC:-Triangularize([p], R);
        
        # Convert lrc to a constructible set
        lrs := map(RC_CST:-RegularSystem, lrc, R);
        return RC_CST:-ConstructibleSet(lrs, R);
    else
        return RC_CST:-GeneralConstruct([p], [], R);
    end if;
    
end proc:


# ------------------------------------------------------------------------------
# init_with_var
#
# Compute a triangular decomposition of p such that p vanishes at every
# point in the triangular decomposition.
#
# INPUT
#   p ...... Polynomial
#    v ...... Variable
#   R ...... Polynomial Ring
#    opts ... A table containing the options (see VanishingConstructibleSet 
#             header)
#
# OUTPUT
#    Same as VanishingConstructibleSet
# ------------------------------------------------------------------------------
init_with_var := proc(p::polynom, v::name, R::TRDring, opts::table, $) :: TRDcs;
    
    local c::{set(polynom),list(polynom)},
          i::posint,
          lrc::list(TRDrc),
          lrs::TRDlrs;
    
    # v must be the greatest variable of R
    if not isGreatestVariable(v, R) then
        error "v must be the greatest variable of R";
    end if;
    
    # Case where p doesn't contain the variable
    if not v in indets(p) then
        return init_no_var(p, R, opts);
    end if;
    
    # Case where p is constant
    if RC:-TRDis_constant(p, R) then
        if p = 0 then
            return RC:-TRDconstructible_set_whole_space(R);
        else
            return RC:-TRDempty_constructible_set();
        end if;
    end if;
    
    # Check for any constant coefficients
    c := {coeffs(p, v)};
    for i to nops(c) do
        if RC:-TRDis_constant(c[i], R) then
            return RC:-TRDempty_constructible_set();
        end if;
    end do;
    
    # Compute square free factorization of coefficients of c
    map[inplace](RC:-TRDGcdFreeFactorization, c, R);
    
    c := convert(c, list);
    
    if opts['lazy'] then
        lrc := RC:-Triangularize(c, R);
        
        # Convert lrc to a constructible set
        lrs := map(RC_CST:-RegularSystem, lrc, R);
        
        return RC_CST:-ConstructibleSet(lrs, R);
    else
        return RC_CST:-GeneralConstruct(c, [], R);
    end if;
    
end proc:

end module:
