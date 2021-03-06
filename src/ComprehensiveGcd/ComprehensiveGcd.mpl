# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ComprehensiveGcd.mpl                                                    #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 24/2017                                                #
#                                                                         #
# Compute the gcd of two parametric univariate polynomials in the sense   #
# of Lazard. Constraints on parameter values can be provided via a        #
# constructible set, regular system or lists of polynomial equality and   #
# inequation constraints.                                                 #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   ComprehensiveGcd(p1, p2, v, R, options)                               #
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
#       [g_i, cof_p1_i, cof_p2_i, rs_i] ... 'outputType' either 'RS' or   #
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
        cleanRS,
        cleanRS_cofactors,
        cleanCS,
        compute_cofactors_rs_list,
        compute_cofactors_rs,
        pseudo_cofactor;

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
    if nargs < 4 then
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
    
    elif RC:-TRDis_polynomial_ring(args[4]) then
        # ComprehensiveGcd(p1, p2, v, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveGcd called as ComprehensiveGcd(p1, p2, v, R, options)");
        
        p1 := args[1];
        p2 := args[2];
        v  := args[3];
        R  := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_F_H(p1, p2, v, [], [], R, opts);
        
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
#       cofactors                                                         #
#       output_CS                                                         #
#       output_RS                                                         #
#    See ComprehensiveGcd header for specifications.                      #
# ----------------------------------------------------------------------- #
processOptions := proc(opts_in::set(equation), $) :: table;

    local opts::table(),
          opt::equation,
          tab_opts,
          opt_name,
          opt_value;

    # Default values
    opts['output_CS'] := true;
    opts['output_RS'] := false;
    opts['cofactors']  := false;
    
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
            elif opt_name = ('cofactors') then
                opt_value := rhs(opt);
                if not type(opts['cofactors'], 'truefalse') then
                    error "incorrect option value: %1", opt;
                end if;
                opts['cofactors'] := opt_value;
            else
                error "unknown option";
            end if;
            tab_opts[opt_name] := true;
        else
            error "incorrect option format";
        end if;
    end do;
    
    # Check option compatibility
    if opts['output_CS'] and opts['cofactors'] then
        error "Output options ConstructibleSet or CS and cofactors=true are not compatible.";
    end if;
    
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
    checkInput(p1, p2, v, R);

    # All elements of F must be polynomials in R
    for i to nops(F) do
        if not RC:-TRDis_poly(F[i], R) then
            error "Invalid polynomial in F";
        end if;
        if v in indets(H[i]) then
            error "Input polynomial equation list F should not contain conditions on %1", v;
        end if;
    end do;

    # All elements of H must be polynomials in R
    for i to nops(H) do
        if not RC:-TRDis_poly(H[i], R) then
            error "Invalid polynomial in H";
        end if;
        if v in indets(H[i]) then
            error "Input polynomial inequation list H should not contain conditions on %1", v;
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
    checkInput(p1, p2, v, R);

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
    checkInput(p1, p2, v, R);

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
# ----------------------------------------------------------------------- #
checkInput := proc(p1::depends(polyInRing(R)), p2::depends(polyInRing(R)), v::name, R::TRDring, $)

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

    # Output options
    if opts['output_RS'] then
        result := convertListWithCSToListWithRS(result, 2, R);
        if opts['cofactors'] then
            result := compute_cofactors_rs_list(p1, p2, result, v, R);
            result := cleanRS_cofactors(result, v, R);
        else
            result := cleanRS(result, v, R);
        end if;
    else
        result := cleanCS(result, R);
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


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# External Files
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveGcd/comprehensive_gcd_src.mpl>
$include <src/ComprehensiveGcd/compute_cofactors_rs.mpl>
$include <src/ComprehensiveGcd/pseudo_cofactor.mpl>
$include <src/ComprehensiveGcd/cleanRS.mpl>
$include <src/ComprehensiveGcd/cleanRS_cofactors.mpl>
$include <src/ComprehensiveGcd/cleanCS.mpl>

end module;