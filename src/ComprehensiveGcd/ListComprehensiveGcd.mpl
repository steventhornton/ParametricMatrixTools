# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ListComprehensiveGcd.mpl                                                #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 24/2017                                                #
#                                                                         #
# Compute the gcd of a list of parametric univariate polynomials in the   #
# sense of Lazard. Constraints on parameter values can be provided via a  #
# constructible set, regular system or lists of polynomial equality and   #
# inequation constraints.                                                 #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   ListComprehensiveGcd(lp, v, R, options)                               #
#   ListComprehensiveGcd(lp, v, rs, R, options)                           #
#   ListComprehensiveGcd(lp, v, cs, R, options)                           #
#   ListComprehensiveGcd(lp, v, F, R, options)                            #
#   ListComprehensiveGcd(lp, v, F, H, R, options)                         #
#                                                                         #
# INPUT                                                                   #
#   lp ... List of polynomials                                            #
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
# OUTPUT                                                                  #
#   A sequence gcdList, cs_zero where gcdList is a list with elements of  #
#   the form:                                                             #
#       [g_i, rs_i] ....................... 'outputType' either 'RS' or   #
#                                           'RegularSystem' and           #
#                                           'cofactors' = false.          #
#       [g_i, cs_i] ....................... 'outputType' either 'CS' or   #
#                                           'ConstructibleSet' and        #
#                                           'cofactors' = false.          #
#   Where g_i is the gcd of the polynomials in lp for all values in the   #
#   zero set of cs_i or rs_i.                                             #
#   cs_zero is the constructible set where both 2 or more polynomials in  #
#   lp vanish for all values in its zero set. The set                     #
#   {cs_1, cs_2, ..., cs_zero} forms a partition of the input             #
#   constructible set.                                                    #
#                                                                         #
# ASSUMPTIONS                                                             #
#   v is the largest variable in R that appears in any polynomial in lp.  #
#   cs must only contain polynomials in variables strictly less than v.   #
#
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
ListComprehensiveGcd := module()

    export ModuleApply;

    local
        init,
        init_F_H,
        init_rs,
        init_cs,
        processOptions,
        checkInput,
        implementation,
        list_comprehensive_gcd_src,
        convertToRS,
        cleanRS,
        cleanCS;

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
#   Same as ListComprehensiveGcd                                          #
# ----------------------------------------------------------------------- #
init := proc()

    local lp, v, F, H, R, rs, cs, opts;

    # Check the number of arguments
    if nargs < 4 then
        error "Insufficient number of arguments";
    elif nargs > 6 then
        error "Too many arguments";
    end if;

    if type(args[3], 'list') and type(args[4], 'list') then
        # ListComprehensiveGcd(p1, p2, v, F, H, R, options)
        userinfo(2, 'ParametricMatrixTools', "ListComprehensiveGcd called as ListComprehensiveGcd(lp, v, F, H, R, options)");

        if nargs = 4 then
            error "Expected a fifth argument of a polynomial ring";
        end if;

        lp := args[1];
        v  := args[2];
        F  := args[3];
        H  := args[4];
        R  := args[5];

        opts := processOptions({args[6..-1]});

        return init_F_H(lp, v, F, H, R, opts);

    elif type(args[3], 'list') then
        # ListComprehensiveGcd(p1, p2, v, F, R, options)
        userinfo(2, 'ParametricMatrixTools', "ListComprehensiveGcd called as ListComprehensiveGcd(lp, v, F, R, options)");

        lp := args[1];
        v  := args[2];
        F  := args[3];
        R  := args[4];

        opts := processOptions({args[5..-1]});

        return init_F_H(lp, v, F, [], R, opts);

    elif RC:-TRDis_regular_system(args[3]) then
        # ListComprehensiveGcd(lp, v, rs, R, options)
        userinfo(2, 'ParametricMatrixTools', "ListComprehensiveGcd called as ListComprehensiveGcd(lp, v, rs, R, options)");

        lp := args[1];
        v  := args[2];
        rs := args[3];
        R  := args[4];

        opts := processOptions({args[5..-1]});

        return init_rs(lp, v, rs, R, opts);


    elif RC:-TRDis_constructible_set(args[3]) then
        # ListComprehensiveGcd(lp, v, cs, R, options)
        userinfo(2, 'ParametricMatrixTools', "ListComprehensiveGcd called as ListComprehensiveGcd(lp, v, cs, R, options)");

        lp := args[1];
        v  := args[2];
        cs := args[3];
        R  := args[4];

        opts := processOptions({args[5..-1]});

        return init_cs(lp, v, cs, R, opts);
    
    elif RC:-TRDis_polynomial_ring(args[3]) then
        # ListComprehensiveGcd(lp, v, R, options)
        userinfo(2, 'ParametricMatrixTools', "ListComprehensiveGcd called as ListComprehensiveGcd(lp, v, R, options)");
        
        lp := args[1];
        v  := args[2];
        R  := args[3];
        
        opts := processOptions({args[4..-1]});
        
        return init_F_H(lp, v, [], [], R, opts);
        
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
#    See ListComprehensiveGcd header for specifications.                  #
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
#   lp ..... List of polynomials                                          #
#   v ...... Variable                                                     #
#   F ...... List of polynomials over R representing equations            #
#   H ...... List of polynomials over R representing inequations          #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see ListComprehensiveGcd     #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ListComprehensiveGcd                                         #
# ----------------------------------------------------------------------- #
init_F_H := proc(lp::depends(list(polyInRing(R))), v::name, F::list(polynom), H::list(polynom), R::TRDring, opts::table, $)

    local i::posint,
          cs::TRDcs;

    # Check the input for errors
    checkInput(lp, v, R);

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

    return implementation(lp, v, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# init_rs                                                                 #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   lp ..... List of polynomials                                          #
#   v ...... Variable                                                     #
#   rs ..... Regular system                                               #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see ListComprehensiveGcd     #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ListComprehensiveGcd                                         #
# ----------------------------------------------------------------------- #
init_rs := proc(lp::depends(list(polyInRing(R))), v::name, rs::TRDrs, R::TRDring, opts::table, $)

    local cs :: TRDcs;

    # Check the input for errors
    checkInput(lp, v, R);

    # All polynomial equations and inequations in rs should be not contain
    # any variables strictly greater than v as an indeterminant.
    if not isUnder(rs, v, R) then
        error "Input regular system should not contain conditions on %1", v;
    end if;

    # Convert rs to a constructible set
    cs := RC_CST:-ConstructibleSet([rs], R);

    return implementation(lp, v, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# init_cs                                                                 #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   lp ..... List of polynomials                                          #
#   v ...... Variable                                                     #
#   cs ..... Constructible set                                            #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see ListComprehensiveGcd     #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ListComprehensiveGcd                                         #
# ----------------------------------------------------------------------- #
init_cs := proc(lp::depends(list(polyInRing(R))), v::name, cs::TRDcs, R::TRDring, opts::table, $)

    # Check the input for errors
    checkInput(lp, v, R);

    # cs should not contain any condition on v
    if not isUnder(cs, v, R) then
        error "Input constructible set should not contain conditions on %1", v;
    end if;

    return implementation(lp, v, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# checkInput                                                              #
#                                                                         #
# Checks for errors in the input values                                   #
#                                                                         #
# INPUT                                                                   #
#   lp ...... List of polynomials                                         #
#   v ...... Variable                                                     #
#   R ...... Polynomial ring                                              #
# ----------------------------------------------------------------------- #
checkInput := proc(lp::depends(list(polyInRing(R))), v::name, R::TRDring, $)
    
    local p :: polynom;
    
    # No polynomials in lp may contain variables strictly greater than v
    for p in lp do
        if not RC:-TRDis_constant(p, R) then
            if RC:-TRDstrictly_less_var(v, RC:-MainVariable(p, R), R) then
                error "No polynomials in lp may contain any variables stricly greater than v";
            end if;
        end if;
    end do;
    
end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Computes the gcd using the specified method and returns the specified   #
# type. Assume no errors in input values.                                 #
#                                                                         #
# INPUT                                                                   #
#   lp_in ... List of polynomials                                         #
#   v ....... Variable                                                    #
#   cs ...... Constructible set                                           #
#   R ....... Polynomial ring                                             #
#   opts .... A table containing the options (see ListComprehensiveGcd    #
#             header)                                                     #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ListComprehensiveGcd                                         #
# ----------------------------------------------------------------------- #
implementation := proc(lp_in::depends(list(polyInRing(R))), v::name, cs::TRDcs, R::TRDring, opts::table, $)

    local lp :: list(polynom),
          result,
          cs_zero :: TRDcs;
    
    lp := map(expand, lp_in);
    
    # Call the algorithm
    result, cs_zero := list_comprehensive_gcd_src(lp, v, cs, R);

    # Convert to regular systems
    if opts['output_RS'] then
        result := convertToRS(result, R);
        result := cleanRS(result, v, R);
    else
        result := cleanCS(result, R);
        
        # Clean result in only one gcd is found
        if nops(result) = 1 and RC:-TRDis_empty_constructible_set(cs_zero, R) then
            result := [[result[1][1], cs]];
        end if;
    end if;

    return result, cs_zero;

end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# External Files
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveGcd/list_comprehensive_gcd_src.mpl>
$include <src/ComprehensiveGcd/cleanRS.mpl>
$include <src/ComprehensiveGcd/cleanCS.mpl>
$include <src/ComprehensiveGcd/convertToRS.mpl>

end module;