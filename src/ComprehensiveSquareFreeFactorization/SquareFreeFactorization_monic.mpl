# ======================================================================= #
# ======================================================================= #
#                                                                         #
# SquareFreeFactorization_monic.mpl                                       #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Mar. 9/2017                                                 #
#                                                                         #
# Compute the square-free factorization of a parametric, univariate       #
# polynomial that is monic in its main variable. A complete case          #
# discussion forming a partition of the input reular system is returned.  #
# The partition of the input regular system is such that over each branch #
# the square-free factorization of the input polynomial is continuous.    #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   SquareFreeFactorization_monic(p, v, R, options)                       #
#   SquareFreeFactorization_monic(p, v, rs, R, options)                   #
#   SquareFreeFactorization_monic(p, v, cs, R, options)                   #
#   SquareFreeFactorization_monic(p, v, F, R, options)                    #
#   SquareFreeFactorization_monic(p, v, F, H, R, options)                 #
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
#   lcoeff(p, mvar(R)) = 1                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list of with elements of the form                                   #
#       [lp_i, cs_i]                                                      #
#   or                                                                    #
#       [lp_i, rs_i]                                                      #
#   where                                                                 #
#       - cs_i is a constructible set (outputType = 'CS')                 #
#       - rs_i is a regular system    (outputType = 'RS')                 #
#       - lp_i is a list with elements of the form:                       #
#             [p_j, n_j]                                                  #
#    such that p = m_i*product(p_j^n_j) and p_j are the square-free       #
#    factors in the zero set of cs_i or rs_i and m_i is some rational     #
#    function of the paramters.                                           #
#                                                                         #
# EXAMPLE                                                                 #
#   > p := (x+1)^2*(x+a):                                                 #
#   > R := PolynomialRing([x, a]):                                        #
#   > sqrFreeFac := SquareFreeFactorization_monic(p, x, R):               #
#   > Display(sqrFreeFac[1], R)                                           #
#         [[[x+a, 1], [x+1, 2]], a-1 <> 0]                                #
#   > Display(sqrFreeFac[2], R)                                           #
#         [[[x+1, 3]], a-1 = 0]                                           #
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
SquareFreeFactorization_monic := module()
    
    export ModuleApply;
    
    local init,
          init_F_H,
          init_rs,
          init_cs,
          processOptions,
          checkInput,
          implementation,
          convertToCS,
          square_free_factorization_monic;
    
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
#   Same as SquareFreeFactorization_monic                                 #
# ----------------------------------------------------------------------- #
init := proc()

    local p, v, F, H, R, rs, cs, opts;

    # Check the number of arguments
    if nargs < 3 then
        error "Insufficient number of arguments";
    elif nargs > 7 then
        error "Too many arguments";
    end if;

    if type(args[3], 'list') and type(args[4], 'list') then
        # SquareFreeFactorization_monic(p, v, F, H, R, options)
        userinfo(2, 'ParametricMatrixTools', "SquareFreeFactorization_monic called as SquareFreeFactorization_monic(p, v, F, H, R, options)");

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
        # SquareFreeFactorization_monic(p, v, F, R, options)
        userinfo(2, 'ParametricMatrixTools', "SquareFreeFactorization_monic called as SquareFreeFactorization_monic(p, v, F, R, options)");

        p := args[1];
        v := args[2];
        F := args[3];
        R := args[4];

        opts := processOptions({args[5..-1]});

        return init_F_H(p, v, F, [], R, opts);

    elif RC:-TRDis_regular_system(args[3]) then
        # SquareFreeFactorization_monic(p, v, rs, R, options)
        userinfo(2, 'ParametricMatrixTools', "SquareFreeFactorization_monic called as SquareFreeFactorization_monic(p, v, rs, R, options)");

        p  := args[1];
        v  := args[2];
        rs := args[3];
        R  := args[4];

        opts := processOptions({args[5..-1]});

        return init_rs(p, v, rs, R, opts);

    elif RC:-TRDis_constructible_set(args[3]) then
        # SquareFreeFactorization_monic(p, v, cs, R, options)
        userinfo(2, 'ParametricMatrixTools', "SquareFreeFactorization_monic called as SquareFreeFactorization_monic(p, v, cs, R, options)");

        p  := args[1];
        v  := args[2];
        cs := args[3];
        R  := args[4];

        opts := processOptions({args[5..-1]});

        return init_cs(p, v, cs, R, opts);
    
    elif RC:-TRDis_polynomial_ring(args[3]) then
        # SquareFreeFactorization_monic(p, v, R, options)
        userinfo(2, 'ParametricMatrixTools', "SquareFreeFactorization_monic called as SquareFreeFactorization_monic(p, v, R, options)");
        
        p  := args[1];
        v  := args[2];
        R  := args[3];
        
        opts := processOptions({args[4..-1]});
        
        return init_F_H(p, v, [], [], R, opts);
        
    else
        error "Invalid arguments";
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
#       output_CS                                                         #
#       output_RS                                                         #
#    See SquareFreeFactorization_monic header for specifications.         #
# ----------------------------------------------------------------------- #
processOptions := proc(opts_in::set(equation), $) :: table;

    local opts :: table(),
          opt :: equation,
          tab_opts,
          opt_name,
          opt_value;

    # Default values
    opts['output_CS'] := true;
    opts['output_RS'] := false;
    
    tab_opts := table();
    
    # Process each option
    for opt in opts_in do
        if type(opt, 'equation') then
            
            opt_name := lhs(opt);
            
            if assigned(tab_opts[opt_name]) then
                error "duplicate option";
            end if;
            
            if opt_name = ('outputType') then
                opt_value := rhs(opt);
                if not member(opt_value, ['CS', 'RS', 'ConstructibleSet', 'RegularSystem']) then
                    error "incorrect option value: %1", opt;
                end if;
                if member(opt_value, ['CS', 'ConstructibleSet']) then
                    opts['output_CS'] := true;
                    opts['output_RS'] := false;
                elif member(opt_value, ['RS', 'RegularSystem']) then
                    opts['output_CS'] := false;
                    opts['output_RS'] := true;
                else
                    error "incorrect option value: %1", opt;
                end if;
            else
                error "unknown option";
            end if;
            tab_opts[opt_name] := true;
        else
            error "incorrect option format";
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
#            SquareFreeFactorization_monic header)                        #
#                                                                         #
# OUTPUT                                                                  #
#    Same as SquareFreeFactorization_monic                                #
# ----------------------------------------------------------------------- #
init_F_H := proc(p::depends(polyInRing(R)), v::name, F::list(polynom), H::list(polynom), R::TRDring, opts::table, $)

    local i :: posint,
          cs :: TRDcs;

    # Check the input for errors
    checkInput(p, v, R);

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
#            SquareFreeFactorization_monic header)                        #
#                                                                         #
# OUTPUT                                                                  #
#    Same as SquareFreeFactorization_monic                                #
# ----------------------------------------------------------------------- #
init_rs := proc(p::depends(polyInRing(R)), v::name, rs::TRDrs, R::TRDring, opts::table, $)

    local cs :: TRDcs;

    # Check the input for errors
    checkInput(p, v, R);

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
#            SquareFreeFactorization_monic header)                        #
#                                                                         #
# OUTPUT                                                                  #
#    Same as SquareFreeFactorization_monic                                #
# ----------------------------------------------------------------------- #
init_cs := proc(p::depends(polyInRing(R)), v::name, cs::TRDcs, R::TRDring, opts::table, $)

    # Check the input for errors
    checkInput(p, v, R);

    # cs should not contain any condition on v
    if not isUnder(cs, v, R) then
        error "Input constructible set should not contain conditions on %1", v;
    end if;

    return implementation(p, v, cs, R, opts);

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
# ----------------------------------------------------------------------- #
checkInput := proc(p::depends(polyInRing(R)), v::name, R::TRDring, $)

    # p1 and p2 must not contain any variables strictly greater than v
    if not RC:-TRDis_constant(p, R) then
        if RC:-TRDstrictly_less_var(v, RC:-MainVariable(p, R), R) then
            error "p must not contain any variables stricly greater than v";
        end if;
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
#   opts .... Table of options, see the                                   #
#             SquareFreeFactorization_monic header for details.           #
#                                                                         #
# OUTPUT                                                                  #
#    Same as SquareFreeFactorization_monic                                #
# ----------------------------------------------------------------------- #
implementation := proc(p_in::depends(polyInRing(R)), v::name, cs::TRDcs, R::TRDring, opts, $)

    local p :: polynom,
          result,
          lrs :: TRDlrs,
          rs :: TRDrs;

    p := expand(p_in);
    
    result := [];
    
    # Call the algorithm
    lrs := RC_CST:-RepresentingRegularSystems(cs, R);
    for rs in lrs do
        result := [op(result), op(square_free_factorization_monic(p, v, rs, R))];
    end do;
    
    if opts['output_CS'] then
        return convertToCS(result, R);
    else
        return result;
    end if;

end proc;


# ----------------------------------------------------------------------- #
# convertToCS                                                             #
#                                                                         #
# Convert the output structure from regular systems to constructible sets.#
#                                                                         #
# INPUT                                                                   #
#   result ...
#   R ........
#                                                                         #
# OUTPUT                                                                  #
#    Same as SquareFreeFactorization_monic                                #
# ----------------------------------------------------------------------- #
convertToCS := proc(result, R, $)
    
    local out, item, cs::TRDcs;
    
    out := [];
    for item in result do
        cs := RC_CST:-ConstructibleSet([item[2]], R);
        out := [op(out), [item[1], cs]];
    end do;
    
    return out;
    
end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# External Files
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveSquareFreeFactorization/square_free_factorization_monic.mpl>

end module;