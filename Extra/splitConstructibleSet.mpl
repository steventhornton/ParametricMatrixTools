# ======================================================================= #
# ======================================================================= #
#
# splitConstructibleSet.mpl
#
# AUTHOR .... Steven E. Thornton
#                Under the supervision of 
#                Robert M. Corless & Marc Moreno Maza
# EMAIL ..... sthornt7@uwo.ca
# UPDATED ... Jan. 18/2016
#
# Split a constructible into two constructible sets such that the polynomial 
# vanishes everywhere on the first constructible set and nowhere on the second.
#
# CALLING SEQUENCE
# 	splitConstructibleSet(p, cs, R)
#   splitConstructibleSet(p, v, cs, R)
#
# INPUT:
#   p .... Polynomial
#	v .... Variable
#   cs ... Constructible set
#   R .... Polynomial ring
#
# OUTPUT:
#	If three arguments are provided:
#		Two constructible sets are returned such that p vanishes at all points 
#		in the first constructible set and at no points exist in the second. The 
#		two constructible sets form a partition of the input constructible set.
#
# 	If four arguments are provided:
#		1. v must be the greatest variable of R
#		2. All polynomials of cs MUST not contain v
#		3. Output is the same as the three agument case except no conditions on
#		   v appear in the output constructible sets.
# ======================================================================= #
# ======================================================================= #
splitConstructibleSet := module()

    export ModuleApply;

    local
        init,
        init_no_var,
        init_with_var,
        implementation_no_var,
        implementation_with_var;

    ModuleApply := proc() :: TRDcs, TRDcs;
        return init(args)
    end proc:


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ------------------------------------------------------------------------------
# init
#
# Calls the appropriate method based on the number of input arguments.
#
# INPUT/OUTPUT
#	Same as splitConstructibleSet
# ------------------------------------------------------------------------------
init := proc() :: TRDcs, TRDcs;

    if nargs = 3 then
        return init_no_var(args);
    elif nargs = 4 then
        return init_with_var(args);
    else
        error "Incorrect number of arguments";
    end if;

end proc:


# ------------------------------------------------------------------------------
# init_no_var
#
# Checks the types of the input and calls the implementation if all input values
# pass checks.
#
# INPUT
#   p_in ... Polynomial
#   cs ..... Constructible set
#   R ...... Polynomial Ring
#
# OUTPUT
#	Same as splitConstructibleSet
# ------------------------------------------------------------------------------
init_no_var := proc(p_in::polynom, cs::TRDcs, R::TRDring, $) :: TRDcs, TRDcs;

    local p::polynom;

    # Ensure p is a polynomial in R
    if not RC:-TRDis_poly(p_in, R) then
        error "invalid polynomial: %1", p_in;
    else
        p := RC:-TRDmodularize_coefficients(p_in, R);
    end if;

    return implementation_no_var(p, cs, R);

end proc:


# ------------------------------------------------------------------------------
# init_with_var
#
# Checks the types of the input and calls the implementation if all input values
# pass checks.
#
# INPUT
#   p_in ... Polynomial
#	v ...... Variable
#   cs ..... Constructible set
#   R ...... Polynomial ring
#
# OUTPUT
#	Same as splitConstructibleSet
# ------------------------------------------------------------------------------
init_with_var := proc(p_in::polynom, v::name, cs::TRDcs, R::TRDring, $) :: TRDcs, TRDcs;

    local p::polynom;

    # Ensure v is the greatest variable of R
    if not isGreatestVariable(v, R) then
        error "v must be the greatest variable of R";
    end if;

    # Ensure p is a polynomial in R
    if not RC:-TRDis_poly(p_in, R) then
        error "invalid polynomial: %1", p_in;
    else
        p := RC:-TRDmodularize_coefficients(p_in, R);
    end if;

    # cs should not contain any condition on v
    if not isUnder(cs, v, R) then
        error "Input constructible set should not contain conditions on %1", v;
    end if;

    return implementation_with_var(p, v, cs, R);

end proc:


# ------------------------------------------------------------------------------
# implementation_no_var
#
# Split a constructible into two constructible sets such that the polynomial 
# vanishes everywhere on the first constructible set and nowhere on the second.
#
# INPUT
#   p .... Polynomial
#   cs ... Constructible set
#   R .... Polynomial ring
#
# OUTPUT
#	Same as splitConstructibleSet
# ------------------------------------------------------------------------------
implementation_no_var := proc(p::polynom, cs::TRDcs, R::TRDring, $) :: TRDcs, TRDcs;

    local cs_zero::TRDcs,
          cs_nonzero::TRDcs;

    # Create constructible set where p = 0
    cs_zero    := VanishingConstructibleSet(p, R, 'lazy'=false);
    cs_zero    := RC_CST:-Intersection(cs, cs_zero, R);

    # Case where p is nonzero
    cs_nonzero := RC_CST:-Difference(cs, cs_zero, R);

    return cs_zero, cs_nonzero;

end proc:


# ------------------------------------------------------------------------------
# implementation_with_var
#
# Split a constructible into two constructible sets such that the polynomial 
# vanishes everywhere on the first constructible set and nowhere on the second.
# The input polynomial is treated as a univariate polynomial in a specified 
# variable, no conditions on this variable can be present in the input 
# constructible set. Neither of the output constructible sets will contain any 
# conditions on this variable.
#
# Steps
#	1. Check for easy cases (p doesn't contain v, p is constant, etc.)
#	2. Extract coefficents of p w.r.t. v
#	3. Compute square free factorization of each coefficents, put in a set
#	4. Compute the constructible set cs_z where all coefficients are zero
#	5. Intersect cs_z with the input constructible set
#	6. Compute constructible set cs_nz where all the coefficients are not 
#	   simultaneously zero
#	7. Return cs_z and cs_nz
#
# INPUT:
#   in_p ... Polynomial
#	v ...... Variable
#   cs ..... Constructible set
#   R ...... Polynomial ring
#
# OUTPUT:
#	Same as splitConstructibleSet
# ------------------------------------------------------------------------------
implementation_with_var := proc(p::polynom, v::name, cs::TRDcs, R::TRDring, $) :: TRDcs, TRDcs;

    local cs_zero::TRDcs,
          cs_nonzero::TRDcs;

    # Case where all coefficients are zero
    cs_zero := VanishingConstructibleSet(p, v, R, 'lazy'=false);
    cs_zero := RC_CST:-Intersection(cs_zero, cs, R);

    # Case where all coefficients are not simultaneously zero
    cs_nonzero := RC_CST:-Difference(cs, cs_zero, R);

    return cs_zero, cs_nonzero;

end proc:

end module: