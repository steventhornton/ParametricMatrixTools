# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ComprehensiveSquareFreeFactorization.mpl                                #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 22/2016                                                #
#                                                                         #
# Compute the square-free decomposition of a parametric univariate        #
# polynomial over a constructible set by a modified version of Yun's      #
# algorithm for parametric univariate polynomials.                        #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   ComprehensiveSquareFreeFactorization(p, v, rs, R, options)            #
#   ComprehensiveSquareFreeFactorization(p, v, cs, R, options)            #
#   ComprehensiveSquareFreeFactorization(p, v, F, R, options)             #
#   ComprehensiveSquareFreeFactorization(p, v, F, H, R, options)          #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   v .... Variable                                                       #
#   rs ... Regular system                                                 #
#   cs ... Constructible set                                              #
#   F .... List of polynomials over R representing equations              #
#   H .... List of polynomials over R representing inequations            #
#   R .... Polynomial ring                                                #
#                                                                         #
# OPTIONS                                                                 #
#   outputType ... 'ConstructibleSet' or 'CS' (default):                  #
#                      - Output will contain constructible sets           #
#                  'RegularSystem' or 'RS':                               #
#                      - Output will contain regular systems              #
#                                                                         #
# ASSUMPTIONS                                                             #
#   degree(p, mvar(R)) > 0                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list of lists of the form                                           #
#       [lp_i, cs_i]                                                      #
#   or                                                                    #
#       [lp_i, rs_i]                                                      #
#   Where                                                                 #
#       - cs_i is a constructible set (outputType = 'CS')                 #
#       - rs_i is a regular system    (outputType = 'RS')                 #
#       - lp_i is a list with elements of the form:                       #
#             [p_j, m]                                                    #
#         such that p = product(p_j^m) and p_j are the square-free        #
#         factors in the zero set of cs_i or rs_i                         #
#                                                                         #
# EXAMPLE                                                                 #
#   > p := (x+1)^2*(x+a):                                                 #
#   > R := PolynomialRing([x, a]):                                        #
#   > cs := GeneralConstruct([], [], R):                                  #
#   > sqrFreeFac := ComprehensiveSquareFreeFactorization(p, x, cs, R):    #
#   > Display(sqrFreeFac[1], R)                                           #
#         [[[x + a, 1], [x + 1, 2]], a-1 <> 0]                            #
#   > Display(sqrFreeFac[2], R)                                           #
#         [[[x + 1, 3]], a-1 = 0]                                         #
#                                                                         #
# REFERENCES                                                              #
#   - Yun, D. Y. (1976, August). On square-free decomposition algorithms. #
#     In Proceedings of the third ACM symposium on Symbolic and algebraic #
#     computation (pp. 26-35). ACM.                                       #
#                                                                         #
# LICENSE                                                                 #
#   This program is free software: you can redistribute it and/or modify  #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation, either version 3 of the License, or     #
#   any later version.                                                    #
#                                                                         #
#   This program is distributed in the hope that it will be useful,       #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
#   You should have received a copy of the GNU General Public License     #
#   along with this program.  If not, see http://www.gnu.org/licenses/.   #
# ======================================================================= #
# ======================================================================= #
ComprehensiveSquareFreeFactorization := module()
    
    export ModuleApply;
    
    local init,
          init_F_H,
          init_rs,
          init_cs,
          processOptions,
          checkInput,
          implementation,
          comprehensive_square_free_factorization_yun;
    
    ModuleApply := proc()
        return init(args);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# ----------------------------------------------------------------------- #
# init                                                                    #
#                                                                         #
# Checks the types of the input and calls the implementation if all input #
# values pass checks.                                                     #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as ComprehensiveSquareFreeFactorization                          #
# ----------------------------------------------------------------------- #
init := proc()

    local p, v, F, H, R, rs, cs, opts;

    # Check the number of arguments
    if nargs < 4 then
        error "Insufficient number of arguments";
    elif nargs > 6 then
        error "Too many arguments";
    end if;

    if type(args[3], 'list') and type(args[2], 'list') then
        # ComprehensiveSquareFreeFactorization(p, v, F, H, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveSquareFreeFactorization called as ComprehensiveSquareFreeFactorization(p, v, F, H, R, options)");

        if nargs = 4 then
            error "Expected a fifth argument of a polynomial ring";
        end if;

        p := args[1];
        v := args[2];
        F := args[3];
        H := args[4];
        R := args[5];

        opts := processOptions({args[6..-1]});

        return init_F_H(p, v, F, H, R, opts);

    elif type(args[3], 'list') then
        # ComprehensiveSquareFreeFactorization(p, v, F, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveSquareFreeFactorization called as ComprehensiveSquareFreeFactorization(p, v, F, R, options)");

        p := args[1];
        v := args[2];
        F := args[3];
        R := args[4];

        opts := processOptions({args[5..-1]});

        return init_F_H(p, v, F, [], R, opts);

    elif RC:-TRDis_regular_system(args[3]) then
        # ComprehensiveSquareFreeFactorization(p, v, rs, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveSquareFreeFactorization called as ComprehensiveSquareFreeFactorization(p, v, rs, R, options)");

        p  := args[1];
        v  := args[2];
        rs := args[3];
        R  := args[4];

        opts := processOptions({args[5..-1]});

        return init_rs(p, v, rs, R, opts);

    elif RC:-TRDis_constructible_set(args[3]) then
        # ComprehensiveSquareFreeFactorization(p, v, cs, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveSquareFreeFactorization called as ComprehensiveSquareFreeFactorization(p, v, cs, R, options)");

        p  := args[1];
        v  := args[2];
        cs := args[3];
        R  := args[4];

        opts := processOptions({args[5..-1]});

        return init_cs(p, v, cs, R, opts);

    else
        error "Expected third argument to be a list of polynomials, regular system or a constructible set";
    end if;

end proc;


# ----------------------------------------------------------------------- #
# processOptions                                                          #
#                                                                         #
# Extracts the options from a set, if an option is missing the default    #
# values is returned.                                                     #
#                                                                         #
# INPUT                                                                   #
#   A set of equations corresponding to the options.                      #
#                                                                         #
# OUTPUT                                                                  #
#    A table with indices                                                 #
#        'outputType'                                                     #
#    See ComprehensiveSquareFreeFactorization header for specifications.  #
# ----------------------------------------------------------------------- #
processOptions := proc(opts_in::set(equation), $) :: table;

    local opts :: table(),
          opt :: equation;

    # Default values
    opts['outputType'] := 'CS';

    # Process each option
    for opt in opts_in do
        if lhs(opt) in {indices(opts, 'nolist')} then
            opts[lhs(opt)] := rhs(opt);
        else
            error "'%1' is not a valid option", lhs(opt);
        end if;
    end do;

    return opts;

end proc;


# ----------------------------------------------------------------------- #
# init_F_H                                                                #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   p ...... Polynomial                                                   #
#   v ...... Variable                                                     #
#   F ...... List of polynomials over R representing equations            #
#   H ...... List of polynomials over R representing inequations          #
#   R ...... Polynomial Ring                                              #
#   opts ... A table containing the options (see                          #
#            ComprehensiveSquareFreeFactorization header)                 #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveSquareFreeFactorization                         #
# ----------------------------------------------------------------------- #
init_F_H := proc(p::depends(polyInRing(R)), v::name, F::list(polynom), H::list(polynom), R::TRDring, opts::table, $)

    local i :: posint,
          cs :: TRDcs;

    # Check the input for errors
    checkInput(p, v, R, opts);

    # All elements of F must be polynomials in R
    for i to nops(F) do
        if not RC:-TRDis_poly(F[i], R) then
            error "Invalid polynomial in F";
        end if;
        if v in indets(F[i]) then
            error "Input polynomial equation list F should not contain conditions on %1", v;
        end if;
    end do;

    # All elements of H must be polynomials in R
    for i to nops(H) do
        if not RC:-TRDis_poly(H[i], R) then
            error "Invalid polynomial in H";
        end if;
        if v in indets(F[i]) then
            error "Input polynomial inequation list H should not contain conditions on %1", v;
        end if;
    end do;

    # Convert F and H to a constructible set
    cs := RC_CST:-GeneralConstruct(F, H, R);

    return implementation(p, v, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# init_rs                                                                 #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   p ...... Polynomial                                                   #
#   v ...... Variable                                                     #
#   rs ..... Regular system                                               #
#   R ...... Polynomial Ring                                              #
#   opts ... A table containing the options (see                          #
#            ComprehensiveSquareFreeFactorization header)                 #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveSquareFreeFactorization                         #
# ----------------------------------------------------------------------- #
init_rs := proc(p::depends(polyInRing(R)), v::name, rs::TRDrs, R::TRDring, opts::table, $)

    local cs :: TRDcs;

    # Check the input for errors
    checkInput(p, v, R, opts);

    # All polynomial equations and inequations in rs should be not contain
    # any variables strictly greater than v as an indeterminant.
    if not isUnder(rs, v, R) then
        error "Input regular system should not contain conditions on %1", v;
    end if;

    # Convert rs to a constructible set
    cs := RC_CST:-ConstructibleSet([rs], R);

    return implementation(p, v, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# init_cs                                                                 #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   p ...... Polynomial                                                   #
#   v ...... Variable                                                     #
#   cs ..... Constructible set                                            #
#   R ...... Polynomial Ring                                              #
#   opts ... A table containing the options (see                          #
#            ComprehensiveSquareFreeFactorization header)                 #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveSquareFreeFactorization                         #
# ----------------------------------------------------------------------- #
init_cs := proc(p::depends(polyInRing(R)), v::name, cs::TRDcs, R::TRDring, opts::table, $)

    # Check the input for errors
    checkInput(p, v, R, opts);

    # cs should not contain any condition on v
    if not isUnder(cs, v, R) then
        error "Input constructible set should not contain conditions on %1", v;
    end if;

    return implementation(p, v, cs, R);

end proc;


# ----------------------------------------------------------------------- #
# checkInput                                                              #
#                                                                         #
# Checks for errors in the input values.                                  #
#                                                                         #
# INPUT                                                                   #
#   p ...... Polynomial                                                   #
#   v ...... Variable                                                     #
#   R ...... Polynomial Ring                                              #
#   opts ... A table containing the options (see                          #
#            ComprehensiveSquareFreeFactorization header)                 #
# ----------------------------------------------------------------------- #
checkInput := proc(p::depends(polyInRing(R)), v::name, R::TRDring, opts::table, $)

    # p1 and p2 must not contain any variables strictly greater than v
    if not RC:-TRDis_constant(p, R) then
        if RC:-TRDstrictly_less_var(v, RC:-MainVariable(p, R), R) then
            error "p must not contain any variables stricly greater than v";
        end if;
    end if;

    # outputType option must be either 'RegularSystem', 'RS', 
    # 'ConstructibleSet', or 'CS'
    if not opts['outputType'] in {'RegularSystem', 'RS', 'ConstructibleSet', 'CS'} then
        error "outputType option must be either RegularSystem, RS, ConstructibleSet or CS";
    end if;
    if opts['outputType'] = 'RegularSystem' then
        opts['outputType'] := 'RS';
    end if;
    if opts['outputType'] = 'ConstructibleSet' then
        opts['outputType'] := 'CS';
    end if;

end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Computes the gcd using the specified method and returns the specified   #
# type. Assume no errors in input values.                                 #
#                                                                         #
# INPUT                                                                   #
#   p ....... Polynomial                                                  #
#   v ....... Variable                                                    #
#   cs ...... Constructible set                                           #
#   R ....... Polynomial Ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveSquareFreeFactorization                         #
# ----------------------------------------------------------------------- #
implementation := proc(p_in::depends(polyInRing(R)), v::name, cs::TRDcs, R::TRDring, $)

    local p :: polynom,
          result;

    p := expand(p_in);

    # Call the algorithm
    result := comprehensive_square_free_factorization_yun(p, v, cs, R);
    
    return result;

end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# External Files
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveSquareFreeFactorization/comprehensive_square_free_factorization_yun.mpl>

end module;