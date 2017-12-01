# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ComprehensiveMinimalPolynomial.mpl                                      #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 24/2017                                                #
#                                                                         #
# Computes a complete case discussion for the minimal polynomial of a     #
# matrix where the entries are multivariate polynomials whose             #
# indeterminants are regarded as parameters. Computation is done over     #
# polynomial equality and inequation constraints on the parameters.       #
#                                                                         #
# CALLING SEQUENCE                                                        #
#    ComprehensiveMinimalPolynomial(A, v, R, options)                     #
#    ComprehensiveMinimalPolynomial(A, v, rs, R, options)                 #
#    ComprehensiveMinimalPolynomial(A, v, cs, R, options)                 #
#    ComprehensiveMinimalPolynomial(A, v, F, R, options)                  #
#    ComprehensiveMinimalPolynomial(A, v, F, H, R, options)               #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   v .... Variable for the returned minimal polynomial                   #
#   rs ... Regular system                                                 #
#   cs ... Constructible set                                              #
#   F .... List of polynomials over R representing equations              #
#   H .... List of polynomials over R representing inequations            #
#   R .... Polynomial ring                                                #
#                                                                         #
# OPTIONS                                                                 #
#   outputType ....... Default = CS                                       #
#                      ConstructibleSet or CS:                            #
#                         - Output will contain constructible sets        #
#                      RegularSystem or RS:                               #
#                           - Output will contain regular systems         #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements in one of the following forms:                   #
#       [pmin, rs] .......... 'outputType' is 'RegularSystem' or 'RS'     #
#       [pmin, cs] .......... 'outputType' is 'ConstructibleSet' or 'CS'  #
#   Where pmin is the minimal polynomial of A for all parameter values    #
#   that belonging to the solution set of cs or rs. Together, all the     #
#   constructible sets or regular systems in the output form a partition  #
#   of the input constructible set or regular system.                     #
#                                                                         #
# EXAMPLE                                                                 #
#                                                                         #
# REFERENCES                                                              #
#                                                                         #
# ======================================================================= #
# ======================================================================= #
ComprehensiveMinimalPolynomial := module()
    
    export ModuleApply;
    
    local
        init,
        init_F_H,
        init_rs,
        init_cs,
        
        processOptions,
        checkInput,
        
        implementation,
        
        # ALGORITHMS
        comprehensive_minimal_polynomial_snf;
    
    ModuleApply := proc()
        return init(args);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ----------------------------------------------------------------------- #
# init                                                                    #
#                                                                         #
# Checks the types of the input and calls the appropriate implementation  #
# if all input values pass checks.                                        #
#                                                                         #
# INPUT/OUTPUT                                                            #
#    Same as ComprehensiveMinimalPolynomial                               #
# ----------------------------------------------------------------------- #
init := proc()

    local A, F, H, R, v, opts, cs, rs;

    # Check the number of arguments
    if nargs < 3 then
        error "Insufficient number of arguments";
    elif nargs > 6 then
        error "To many arguments";
    end if;

    if type(args[3], 'list') and type(args[4], 'list') then
        # MinimalPolynomial(A, v, F, H, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveMinimalPolynomial called as ComprehensiveMinimalPolynomial(A, v, F, H, R, options)");

        if nargs = 4 then
            error "Expected a fifth argument of a polynomial ring";
        end if;

        A := args[1];
        v := args[2];
        F := args[3];
        H := args[4];
        R := args[5];
        
        opts := processOptions({args[6..-1]});
        
        return init_F_H(A, v, F, H, R, opts);

    elif type(args[3], 'list') then
        # MinimalPolynomial(A, v, F, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveMinimalPolynomial called as ComprehensiveMinimalPolynomial(A, v, F, R, options)");

        A := args[1];
        v := args[2];
        F := args[3];
        R := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_F_H(A, v, F, [], R, opts);

    elif RC:-TRDis_regular_system(args[3]) then
        # MinimalPolynomial(A, v, rs, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveMinimalPolynomial called as ComprehensiveMinimalPolynomial(A, v, rs, R, options)");

        A  := args[1];
        v  := args[2];
        rs := args[3];
        R  := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_rs(A, v, rs, R, opts);

    elif RC:-TRDis_constructible_set(args[3]) then
        # MinimalPolynomial(A, v, cs, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveMinimalPolynomial called as ComprehensiveMinimalPolynomial(A, v, cs, R, options)");

        A  := args[1];
        v  := args[2];
        cs := args[3];
        R  := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_cs(A, v, cs, R, opts);
    
    elif RC:-TRDis_polynomial_ring(args[3]) then
        # MinimalPolynomial(A, v, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveMinimalPolynomial called as ComprehensiveMinimalPolynomial(A, v, R, options)");
        
        A := args[1];
        v := args[2];
        R := args[3];
        
        opts := processOptions({args[4..-1]});
        
        return init_F_H(A, v, [], [], R, opts);
        
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
#   A table with indices                                                  #
#       output_CS                                                         #
#       output_RS                                                         #
#   See ComprehensiveMinimalPolynomial header for specifications.         #
# ----------------------------------------------------------------------- #
processOptions := proc(opts_in::set(equation), $) :: table;

    local opts::table(),
          tab_opts,
          opt::equation,
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
#   A ...... Matrix of multivariate polynomials                           #
#   v ...... Variable                                                     #
#   F ...... List of polynomials representing equations                   #
#   H ...... List of polynomials representing inequations                 #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see MinimalPolynomial        #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveMinimalPolynomial                               #
# ----------------------------------------------------------------------- #
init_F_H := proc(A::Matrix, v::name, F::list(polynom), H::list(polynom), R::TRDring, opts::table, $)
    
    local p :: polynom,
          cs :: TRDcs;
    
    # Check the input for errors
    checkInput(A, v, R);
    
    # All elements of F must be polynomials in R
    for p in F do
        if not RC:-TRDis_poly(p, R) then
            error "Invalid polynomial in F";
        end if;
    end do;
    
    # All elements of H must be polynomials in R
    for p in H do
        if not RC:-TRDis_poly(p, R) then
            error "Invalid polynomial in H";
        end if;
    end do;
    
    # Convert F and H to a constructible set
    cs := RC_CST:-GeneralConstruct(F, H, R);
    
    return implementation(A, v, cs, R, opts);
    
end proc;


# ----------------------------------------------------------------------- #
# init_rs                                                                 #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   v ...... Variable                                                     #
#   rs ..... Regular system                                               #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see                          #
#            ComprehensiveMinimalPolynomial header)                       #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveMinimalPolynomial                               #
# ----------------------------------------------------------------------- #
init_rs := proc(A::Matrix, v::name, rs::TRDrs, R::TRDring, opts::table, $)
    
    local cs :: TRDcs;
    
    # Check the input for errors
    checkInput(A, v, R);
    
    # rs must be a regular system
    if not RC:-TRDis_regular_system(rs, R) then
        error "Expected a regular system";
    end if;
    
    # Convert rs to a constructible set
    cs := RC_CST:-ConstructibleSet([rs], R);
    
    return implementation(A, v, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# init_cs                                                                 #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   v ...... Variable                                                     #
#   cs ..... Constructible set                                            #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see                          #
#            ComprehensiveMinimalPolynomial header)                       #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveMinimalPolynomial                               #
# ----------------------------------------------------------------------- #
init_cs := proc(A::Matrix, v::name, cs::TRDcs, R::TRDring, opts, $)

    # Check the input for errors
    checkInput(A, v, R);

    # cs must be a constructible set
    if not RC:-TRDis_constructible_set(cs, R) then
        error "Expected a constructible set";
    end if;
    
    return implementation(A, v, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# checkInput                                                              #
#                                                                         #
# Checks for errors in the input values.                                  #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   v ...... Variable                                                     #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see                          #
#            ComprehensiveMinimalPolynomial header)                       #
# ----------------------------------------------------------------------- #
checkInput := proc(A::Matrix, v::name, R::TRDring, $)
    
    local i :: posint, 
          j :: posint, 
          n :: nonnegint,
          m :: nonnegint;
    
    n, m := LA:-Dimension(A);
    
    if n <> m then
        error "Matrix must be square";
    end if;
    
    # A must be at least a 1x1 matrix
    if n < 1 or m < 1 then
        error "Matrix must not be empty";
    end if;
    
    # All elements of A must be polynomials in R
    # None of the elements of A can be polynomial in v
    for i to n do
        for j to m do
            if not RC:-TRDis_poly(A[i,j], R) then
                error "Expected a matrix of polynomials in R";
            end if;
            if v in indets(A[i,j]) then
                error "Matrix must not contain input variable";
            end if;
        end do;
    end do;
    
    # v should no be a variable in the polynomial ring
    # Note: this will also ensure F, H, cs, and rs do not depend on v
    if v in R['variables'] then
        error "Polynomial ring must not contain input variable";
    end if;
    
end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Computes the minimal polynomial using the specified method and returns  #
# the specified type. Assume no errors in input values.                   #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   v ...... Variable                                                     #
#   cs ..... Constructible set                                            #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see                          #
#            ComprehensiveMinimalPolynomial header)                       #
#                                                                         #
# OUTPUT                                                                  #
#    Same as MinimalPolynomial                                            #
# ----------------------------------------------------------------------- #
implementation := proc(AA::Matrix, v::name, cs::TRDcs, R::TRDring, opts::table, $)
    
    local rs :: TRDrs,
          rc :: TRDrc,
          A :: 'Matrix'(square),
          pmin :: polynom,
          result := [],
          lrsCompute :: TRDlrs := [],
          csCompute :: TRDcs;
    
    # Check for zero or constant matrices for each regular system in cs
    for rs in RC_CST:-RepresentingRegularSystems(cs, R) do
        
        # Map SparsePseudoRemainder to A
        rc := RC_CST:-RepresentingChain(rs, R);
        A := map(RC:-SparsePseudoRemainder, AA, rc, R);
        
        # Check if A is a constant matrix
        if isConstantMatrix(A, R) then
            pmin := LS:-MinimalPolynomial(A, v);
            result := [op(result), [pmin, rs]];
            
        # Check if A is a zero matrix over rs
        elif isZeroMatrixOverRS(A, rs, R) then
            result := [op(result), [v, rs]];
        else
            lrsCompute := [op(lrsCompute), rs];
        end if;    
        
    end do;
    
    # Convert all regular systems in result list to constructible sets
    result := map(input -> [op(input[1..-2]), RC_CST:-ConstructibleSet([input[-1]], R)], result);
    
    # Convert lrsCompute to a constructible set
    csCompute := RC_CST:-ConstructibleSet(lrsCompute, R);
    
    # Call the algorithm
    if not RC_CST:-IsEmpty(csCompute, R) then
        result := [op(result), op(comprehensive_minimal_polynomial_snf(AA, v, csCompute, R))];
    end if;
    
    # Output options
    if opts['output_RS'] then
        result := convertListWithCSToListWithRS(result, 2, R);
    end if;
    
    return result;
    
end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# EXTERNAL FILES
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveMinimalPolynomial/comprehensive_minimal_polynomial_snf.mpl>

end module;