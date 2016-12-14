# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ComprehensiveGcd.mpl                                                    #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 14/2016                                                #
#                                                                         #
# Compute the gcd of two parametric univariate polynomials in the sense   #
# of Lazard. Constraints on parameter values can be provided via a        #
# constructible set, regular system or lists of polynomial equality and   #
# inequation constraints.                                                 #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   ComprehensiveGcd(p1, p2, v, rs, R, options)                           #
#   ComprehensiveGcd(p1, p2, v, cs, R, options)                           #
#   ComprehensiveGcd(p1, p2, v, F, R, options)                            #
#   ComprehensiveGcd(p1, p2, v, F, H, R, options)                         #
#                                                                         #
# INPUT                                                                   #
#   p1 ... Polynomial                                                     #
#   p2 ... Polynomial                                                     #
#   v .... Variable                                                       #
#   rs ... Regular system                                                 #
#   cs ... Constructible set                                              #
#   F .... List of polynomials over R representing equations              #
#   H .... List of polynomials over R representing inequations            #
#   R .... Polynomial ring                                                #
#                                                                         #
# OPTIONS                                                                 #
#   cofactors .... false (default):                                       #
#                     - No cofactors are computed                         #
#                  true:                                                  #
#                     - Cofactors are computed (see output)               #
#   outputType ... 'ConstructibleSet' or 'CS' (default):                  #
#                      - Output will contain constructible sets           #
#                  'RegularSystem' or 'RS':                               #
#                      - Output will contain regular systems              #
#                                                                         #
# OPTION COMPATIBILITY                                                    #
#   - 'cofactors' = true and 'outputType' = 'ConstructibleSet' or 'CS'    #
#     are incompatible.                                                   #
#                                                                         #
# OUTPUT                                                                  #
#   A sequence gcdList, cs_zero where gcdList is a list with elements of  #
#   the form:                                                             #
#       [g_i, rs_i] ....................... 'outputType' either 'RS' or   #
#                                           'RegularSystem' and           #
#                                           'cofactors' = false.          #
#       [g_i, cs_i] ....................... 'outputType' either 'CS' or   #
#                                           'ConstructibleSet' and        #
#                                           'cofactors' = false.          #
#       [c_i, cof_p1_i, cof_p2_i, rs_i] ... 'outputType' either 'RS' or   #
#                                           'RegularSystem' and           #
#                                           'cofactors' = true.           #
#   Where g_i is the gcd of p1 and p2 for all values in the zero set of   #
#   cs_i or rs_i, and cof_p1_i and cof_p2_i are the cofactors of p1 and   #
#   p2 respectively:                                                      #
#       cof_p1_i = p1/g_i                                                 #
#       cof_p2_i = p2/g_i                                                 #
#   cs_zero is the constructible set where both p1 and p2 vanish for all  #
#   values in its zero set. The set {cs_1, cs_2, ..., cs_zero} forms a    #
#   partition of the input constructible set.                             #
#                                                                         #
# ASSUMPTIONS                                                             #
#   v is the largest variable in R that appears in p1 and p2.             #
#   cs must only contain polynomials in variables strictly less than v    #
#                                                                         #
# EXAMPLE                                                                 #
#   > p1 := (x+1)*(x+a):                                                  #
#   > p2 := (x+2)^2:                                                      #
#   > R := PolynomialRing([x, a]):                                        #
#   > cs := GeneralConstruct([], [], R):                                  #
#   > g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R):                   #
#   > IsEmpty(cs_zero, R);                                                #
#         true                                                            #
#   > Display(g, R)                                                       #
#         [[-a*x-2*a+2*x+4, a-2 <> 0], [x^2+4*x+4, a-2 = 0]]              #
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
ComprehensiveGcd := module()

    export ModuleApply;

    local
        init,
        init_F_H,
        init_rs,
        init_cs,
        processOptions,
        checkInput,
        implementation,
        comprehensive_gcd_src,
        convertToRS,
        cleanRS,
        compute_cofactors_rs_list,
        compute_cofactors_rs,
        pseudo_cofactor,
        listGcd;

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
#   Same as ComprehensiveGcd                                              #
# ----------------------------------------------------------------------- #
init := proc()

    local p1, p2, v, F, H, R, rs, cs, opts;

    # Check the number of arguments
    if nargs < 5 then
        error "Insufficient number of arguments";
    elif nargs > 8 then
        error "Too many arguments";
    end if;

    if type(args[4], 'list') and type(args[5], 'list') then
        # ComprehensiveGcd(p1, p2, v, F, H, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveGcd called as ComprehensiveGcd(p1, p2, v, F, H, R, options)");

        if nargs = 5 then
            error "Expected a sixth argument of a polynomial ring";
        end if;

        p1 := args[1];
        p2 := args[2];
        v  := args[3];
        F  := args[4];
        H  := args[5];
        R  := args[6];

        opts := processOptions({args[7..-1]});

        return init_F_H(p1, p2, v, F, H, R, opts);

    elif type(args[4], 'list') then
        # ComprehensiveGcd(p1, p2, v, F, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveGcd called as ComprehensiveGcd(p1, p2, v, F, R, options)");

        p1 := args[1];
        p2 := args[2];
        v  := args[3];
        F  := args[4];
        R  := args[5];

        opts := processOptions({args[6..-1]});

        return init_F_H(p1, p2, v, F, [], R, opts);

    elif RC:-TRDis_regular_system(args[4]) then
        # ComprehensiveGcd(p1, p2, v, rs, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveGcd called as ComprehensiveGcd(p1, p2, v, rs, R, options)");

        p1 := args[1];
        p2 := args[2];
        v  := args[3];
        rs := args[4];
        R  := args[5];

        opts := processOptions({args[6..-1]});

        return init_rs(p1, p2, v, rs, R, opts);


    elif RC:-TRDis_constructible_set(args[4]) then
        # ComprehensiveGcd(p1, p2, v, cs, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveGcd called as ComprehensiveGcd(p1, p2, v, cs, R, options)");

        p1 := args[1];
        p2 := args[2];
        v  := args[3];
        cs := args[4];
        R  := args[5];

        opts := processOptions({args[6..-1]});

        return init_cs(p1, p2, v, cs, R, opts);

    else
        error "Expected fourth argument to be a list of polynomials, regular system or a constructible set";
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
#        'cofactors'                                                      #
#        'outputType'                                                     #
#    See ComprehensiveGcd header for specifications.                      #
# ----------------------------------------------------------------------- #
processOptions := proc(opts_in::set(equation), $) :: table;

    local opts::table(),
          opt::equation;

    # Default values
    opts['outputType'] := 'CS';
    opts['cofactors']  := false;

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
#   p1 ..... Polynomial                                                   #
#   p2 ..... Polynomial                                                   #
#   v ...... Variable                                                     #
#   F ...... List of polynomials over R representing equations            #
#   H ...... List of polynomials over R representing inequations          #
#   R ...... Polynomial Ring                                              #
#   opts ... A table containing the options (see ComprehensiveGcd header) #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveGcd                                             #
# ----------------------------------------------------------------------- #
init_F_H := proc(p1::polynom, p2::polynom, v::name, F::list(polynom), H::list(polynom), R::TRDring, opts::table, $)

    local i::posint,
          cs::TRDcs;

    # Check the input for errors
    checkInput(p1, p2, v, R, opts);

    # All elements of F must be polynomials in R
    for i to nops(F) do
        if not RC:-TRDis_poly(F[i], R) then
            error "Invalid polynomial in F";
        end if;
    end do;

    # All elements of H must be polynomials in R
    for i to nops(H) do
        if not RC:-TRDis_poly(H[i], R) then
            error "Invalid polynomial in H";
        end if;
    end do;

    # Convert F and H to a constructible set
    cs := RC_CST:-GeneralConstruct(F, H, R);

    return implementation(p1, p2, v, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# init_rs                                                                 #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   p1 ..... Polynomial                                                   #
#   p2 ..... Polynomial                                                   #
#   v ...... Variable                                                     #
#   rs ..... Regular system                                               #
#   R ...... Polynomial Ring                                              #
#   opts ... A table containing the options (see ComprehensiveGcd header) #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveGcd                                             #
# ----------------------------------------------------------------------- #
init_rs := proc(p1::polynom, p2::polynom, v::name, rs::TRDrs, R::TRDring, opts::table, $)

    local cs::TRDcs;

    # Check the input for errors
    checkInput(p1, p2, v, R, opts);

    # All polynomial equations and inequations in rs should be not contain
    # any variables strictly greater than v as an indeterminant.
    if not isUnder(rs, v, R) then
        error "Input regular system should not contain conditions on %1", v;
    end if;

    # Convert rs to a constructible set
    cs := RC_CST:-ConstructibleSet([rs], R);

    return implementation(p1, p2, v, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# init_cs                                                                 #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   p1 ..... Polynomial                                                   #
#   p2 ..... Polynomial                                                   #
#   v ...... Variable                                                     #
#   cs ..... Constructible set                                            #
#   R ...... Polynomial Ring                                              #
#   opts ... A table containing the options (see ComprehensiveGcd header) #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveGcd                                             #
# ----------------------------------------------------------------------- #
init_cs := proc(p1::polynom, p2::polynom, v::name, cs::TRDcs, R::TRDring, opts::table, $)

    # Check the input for errors
    checkInput(p1, p2, v, R, opts);

    # cs should not contain any condition on v
    if not isUnder(cs, v, R) then
        error "Input constructible set should not contain conditions on %1", v;
    end if;

    return implementation(p1, p2, v, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# checkInput                                                              #
#                                                                         #
# Checks for errors in the input values                                   #
#                                                                         #
# INPUT                                                                   #
#   p1 ..... Polynomial                                                   #
#   p2 ..... Polynomial                                                   #
#   v ...... Variable                                                     #
#   R ...... Polynomial Ring                                              #
#   opts ... A table containing the options (see ComprehensiveGcd header) #
# ----------------------------------------------------------------------- #
checkInput := proc(p1::depends(polyInRing(R)), p2::depends(polyInRing(R)), v::name, R::TRDring, opts::table, $)

    # p1 and p2 must not contain any variables strictly greater than v
    if not RC:-TRDis_constant(p1, R) then
        if RC:-TRDstrictly_less_var(v, RC:-MainVariable(p1, R), R) then
        error "p1 must not contain any variables stricly greater than v";
        end if;
    end if;
    if not RC:-TRDis_constant(p2, R) then
        if RC:-TRDstrictly_less_var(v, RC:-MainVariable(p2, R), R) then
            error "p2 must not contain any variables stricly greater than v";
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

    # Check the cofactors option
    if not type(opts['cofactors'], 'truefalse') then
        error "cofactors option must be a boolean values";
    end if;

    # Output options 'cofactors' = true and 'outputType' = 'ConstructibleSet'
    # are not compatible
    if opts['outputType'] ='CS' and opts['cofactors'] then
        error "Output options ConstructibleSet or CS and cofactors=true are not compatible.";
    end if;

end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Computes the gcd using the specified method and returns the specified   #
# type. Assume no errors in input values.                                 #
#                                                                         #
# INPUT                                                                   #
#   p1_in ... Polynomial                                                  #
#   p2_in ... Polynomial                                                  #
#   v ....... Variable                                                    #
#   cs ...... Constructible set                                           #
#   R ....... Polynomial Ring                                             #
#   opts .... A table containing the options (see ComprehensiveGcd header)#
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveGcd                                             #
# ----------------------------------------------------------------------- #
implementation := proc(p1_in::depends(polyInRing(R)), p2_in::depends(polyInRing(R)), v::name, cs::TRDcs, R::TRDring, opts, $)

    local p1 :: polynom,
          p2 :: polynom,
          result,
          cs_zero :: TRDcs;

    p1 := expand(p1_in);
    p2 := expand(p2_in);

    # Call the algorithm
    result, cs_zero := comprehensive_gcd_src(p1, p2, v, cs, R);

    # Convert to regular systems
    if opts['outputType'] = 'RS' then
        result := convertToRS(result, R);
        result := cleanRS(result, v, R);
    end if;

    # Compute the cofactors
    if opts['cofactors'] then
        result := compute_cofactors_rs_list(p1, p2, result, v, R);
    end if;

    return result, cs_zero;

end proc;


# ----------------------------------------------------------------------- #
# compute_cofactors_rs_list                                               #
#                                                                         #
# Given a list with elements of the form                                  #
#   [g, rs]                                                               #
# compute the cofactors of p1 and p2 (g = gcd(p1, p2) mod rs) for each    #
# element in the list.                                                    #
#                                                                         #
# INPUT                                                                   #
#   p1 ....... Polynomial                                                 #
#   p2 ....... Polynomial                                                 #
#   result ... A list with elements of the form                           #
#                  [g, rs]                                                #
#              such that g = gcd(p1, p2) for all points in the zero set   #
#              of rs.                                                     #
#   v ........ Variable                                                   #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [g, cof_p1, cof_p2, rs]                                           #
#   with the same number of elements as the input list, and in the same   #
#   order as the input list. g and rs are unchanged. cof_p1 = p1/g and    #
#   cof_2 = p2/g.                                                         #
# ----------------------------------------------------------------------- #
compute_cofactors_rs_list := proc(p1::depends(polyInRing(R)), p2::depends(polyInRing(R)), result, v::name, R::TRDring, $)

    local g :: polynom,
          rs :: TRDrs,
          pair :: [polynom, TRDrs],
          output :: {[],list([polynom, ratpoly, ratpoly, TRDrs])},
          cofactors :: [ratpoly, ratpoly];

    output := [];

    for pair in result do
        g, rs := op(pair);
        cofactors := compute_cofactors_rs(p1, p2, g, v, rs, R);
        output := [op(output), [g, op(cofactors), rs]];
    end do;

end proc;


# ----------------------------------------------------------------------- #
# cleanRS                                                                 #
#                                                                         #
# Given a list with elements of the form                                  #
#   [g, rs]                                                               #
# clean the polynomial g by dividing by its leading coefficient w.r.t v.  #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [g, rs]                                                #
#              where g is a polynomial and rs is a regular system.        #
#   v ........ Variable                                                   #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [g, rs]                                                           #
#   with the same number of elements as the input list, and in the same   #
#   order as the input list. g has been divided by it's leading           #
#   coefficient.                                                          #
# ----------------------------------------------------------------------- #
cleanRS := proc(result::list([polynom, TRDrs]), v::name, R::TRDring, $)

    local output,
          pair :: [polynom, TRDrs],
          g :: ratpoly,
          rs :: TRDrs,
          rc :: TRDrc, 
          inv,
          zdiv, 
          gFactors,
          gMonic, 
          rc_inv::TRDrc;

    output := [];

    for pair in result do
        g, rs := op(pair);

        rc := RC_CST:-RepresentingChain(rs, R);
        
        g := RC:-SparsePseudoRemainder(g, rc, R); 
        
        g := normal(g/listGcd([coeffs(g, v)]));
        
        gMonic := normal(g/lcoeff(g, v));
        
        if denom(gMonic) <> 1 then
            inv, zdiv := op(RC:-Inverse(lcoeff(g, v), rc, R));
            rc_inv := inv[1][3];
            
            ASSERT(nops(zdiv) = 0, "Must not be any zero divisors");
            ASSERT(nops(inv) = 1, "Inverse must only contain one case.");
            ASSERT(RC_CT:-EqualSaturatedIdeals(rc, rc_inv, R), "Inverse must not lose any cases.");
            
            g := RC:-SparsePseudoRemainder(g*inv[1][1], rc, R);
            if denom(normal(g/lcoeff(g, v))) = 1 then
                g := normal(g/lcoeff(g, v));
            end if;
        end if;
        
        g := RC:-SparsePseudoRemainder(g, rc, R); 
        gFactors := factors(g)[2];
        gFactors := remove((x,v) -> not v in indets(x[1]), gFactors, v);
        g := mul(map(x -> x[1]^x[2], gFactors));
        
        output := [op(output), [g, rs]];
        
    end do;

    return output;

end proc;


# ----------------------------------------------------------------------- #
# convertToRS                                                             #
#                                                                         #
# Given a list with elements of the form                                  #
#   [g, cs]                                                               #
# expand the list by extracting all regular systems from each             #
# constructible set.                                                      #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [g, cs]                                                #
#              where g is a polynomial and rs is a regular system.        #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [g, rs]                                                           #
# ----------------------------------------------------------------------- #
convertToRS := proc(result, R::TRDring, $)
    
    local output:: {[], list([polynom, TRDrs])},
          pair :: [polynom, TRDcs],
          g :: polynom,
          cs :: TRDcs,
          lrs :: TRDlrs;
    
    output := [];
    
    for pair in result do
        g, cs := op(pair);
        lrs := RC_CST:-RepresentingRegularSystems(cs, R);
        output := [op(output), op(zip((x, y) -> [x, y], g, lrs))];
    end do;
    
    return output;
    
end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# External Files
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveGcd/comprehensive_gcd_src.mpl>
$include <src/ComprehensiveGcd/comprehensive_gcd_src.mpl>
$include <src/ComprehensiveGcd/comprehensive_gcd_src.mpl>
$include <src/ComprehensiveGcd/compute_cofactors_rs.mpl>
$include <src/ComprehensiveGcd/pseudo_cofactor.mpl>
$include <src/ComprehensiveGcd/listGcd.mpl>

end module;