# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ComprehensiveGcd.mpl                                                    #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 8/2016                                                 #
#                                                                         #
# Computes the gcd of two parametric univariate polynomials in the sense  #
# of Lazard. That is, two univariate polynomials where the coefficients   #
# are multivariate polynomials in the parameters. Computation is done     #
# modulo a regular system or constructible set.                           #
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
#   lazy ......... true (default):                                        #
#                      - Only compute one branch of the computation       #
#                  false:                                                 #
#                      - Compute all branches of the computation          #
#   outputType ... "ConstructibleSet" or "CS" (default):                  #
#                      - Output will contain constructible sets           #
#                  "RegularSystem" or "RS":                               #
#                      - Output will contain regular systems              #
#                                                                         #
# OPTION COMPATIBILITY                                                    #
#   - 'cofactors' = true and 'outputType' = "ConstructibleSet" or "CS"    #
#     are incompatible.                                                   #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements in one of the following forms:                   #
#       [g_i, rs_i] ....................... 'outputType' either 'RS' or   #
#                                           'RegularSystem' and           #
#                                           'cofactors' = false.          #
#       [g_i, cs_i] ....................... 'outputType' either 'CS' or   #
#                                           'ConstructibleSet' and        #
#                                           'cofactors' = false.          #
#       [c_i, cof_p1_i, cof_p2_i, rs_i] ... 'outputType' either 'RS' or   #
#                                           'RegularSystem' and           #
#   Where g_i is the gcd of p1 and p2 for all parameter values that       #
#   satisfy the equations and inequations of cs_i or rs_i. cof_p1_i and   #
#   cof_p2_i are the cofactors of p1 and p2 respectively.                 #
#       cof_p1_i = p1/g                                                   #
#       cof_p2_i = p2/g                                                   #
#   ** If g=0 then cof_p1 and cof_p2 will be zero.                        #
#   Together, all the constructible sets cs_i or regular systems rs_i     #
#   (for all values of i) form a partition of the input constructible set #
#   or regular system.                                                    #
#                                                                         #
# ASSUMPTIONS                                                             #
#                                                                         #
# EXAMPLE                                                                 #
#                                                                         #
# REFERENCES                                                              #
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
        comprehensive_gcd_src;
    
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
# init := proc() :: {list([polynom, TRDrs]), list([polynom, TRDcs]), list([polynom, ratpoly, ratpoly, TRDrs])};
init := proc()
    
    local p1, p2, v, F, H, R, rs, cs, opts;
    
    # Check the number of arguments
    if nargs < 5 then
        error "Insufficient number of arguments";
    elif nargs > 9 then
        error "To many arguments";
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
#        'lazy'                                                           #
#        'outputType'                                                     #
#    See ComprehensiveGcd header for specifications.                      #
# ----------------------------------------------------------------------- #
processOptions := proc(opts_in::set(equation), $) :: table;

    local opts::table(),
          opt::equation;

    # Default values
    opts['outputType'] := "CS";
    opts['lazy']       := true;
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
init_F_H := proc(p1::polynom, p2::polynom, v::name, F::list(polynom), H::list(polynom), R::TRDring, opts::table, $) :: {list([polynom, TRDrs]), list([polynom, TRDcs]), list([polynom, ratpoly, ratpoly, TRDrs])};

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
# init_rs := proc(p1::polynom, p2::polynom, v::name, rs::TRDrs, R::TRDring, opts::table, $) :: {list([polynom, TRDrs]), list([polynom, TRDcs]), list([polynom, ratpoly, ratpoly, TRDrs])};
init_rs := proc(p1::polynom, p2::polynom, v::name, rs::TRDrs, R::TRDring, opts::table, $)

    local cs::TRDcs;
    
    # Check the input for errors
    checkInput(p1, p2, v, R, opts);
    
    # rs should not contain any condition on v
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
# init_cs := proc(p1::polynom, p2::polynom, v::name, cs::TRDcs, R::TRDring, opts::table, $) :: {list([polynom, TRDrs]), list([polynom, TRDcs]), list([polynom, ratpoly, ratpoly, TRDrs])};
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
checkInput := proc(p1::polynom, p2::polynom, v::name, R::TRDring, opts::table, $)
    
    # v must be the greatest variable of R
    if not isGreatestVariable(v, R) then
        error "v must be the greatest variable of R";
    end if;
    
    # p1 must be a polynomial in R
    if not RC:-TRDis_poly(p1, R) then
        error "invalid polynomial: %1", p1;
    end if;
    
    # p2 must be a polynomial in R
    if not RC:-TRDis_poly(p2, R) then
        error "invalid polynomial: %1", p2;
    end if;
    
    # outputType option must be either "RegularSystem", "RS", "ConstructibleSet" 
    # or "CS"
    if not opts['outputType'] in {"RegularSystem", "RS", "ConstructibleSet", "CS"} then
        error "outputType option must be either RegularSystem, RS, ConstructibleSet or CS";
    end if;
    
    # Check the cofactors option
    if not type(opts['cofactors'], 'truefalse') then
        error "cofactors option must be a boolean values";
    end if;
    
    # Check the lazy option
    if not type(opts['lazy'], 'truefalse') then
        error "lazy option must be a boolean values";
    end if;
    
    # Output options 'cofactors' = true and 'outputType' = "ConstructibleSet"
    # are not compatible
    if (opts['outputType'] in {"ConstructibleSet", "CS"}) and opts['cofactors'] then
        error "Output options ConstructibleSet or CS and cofactors=true are not compatible.";
    end if;
    
end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Computes the Gcd using the specified method and returns the specified   #
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
# implementation := proc(p1::polynom, p2::polynom, v::name, cs::TRDcs, R::TRDring, opts::table, $) :: {list([polynom, TRDrs]), list([polynom, TRDcs]), list([polynom, ratpoly, ratpoly, TRDrs])};
implementation := proc(p1_in::polynom, p2_in::polynom, v::name, cs::TRDcs, R::TRDring, opts::table, $)
    
    local p1, p2;
    
    p1 := expand(p1_in);
    p2 := expand(p2_in);
    
    # Call the algorithm
    return comprehensive_gcd_src(p1, p2, v, cs, R);
    
end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# External Files
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <ComprehensiveGcd/comprehensive_gcd_src.mpl>

end module;