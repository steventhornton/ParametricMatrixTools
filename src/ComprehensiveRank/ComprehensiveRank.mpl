# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ComprehensiveRank.mpl                                                   #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 24/2017                                                #
#                                                                         #
# Computes a complete case discussion of the rank of a matrix where the   #
# entries are multivariate polynomials whose indeterminants are treated   #
# as parameters. Computation is done modulo a regular system or           #
# constructible set.                                                      #
#                                                                         #
# CALLING SEQUENCE                                                        #
#    ComprehensiveRank(A, R, options)                                     #
#    ComprehensiveRank(A, rs, R, options)                                 #
#    ComprehensiveRank(A, cs, R, options)                                 #
#    ComprehensiveRank(A, F, R, options)                                  #
#    ComprehensiveRank(A, F, H, R, options)                               #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
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
#       [r, rs] .......... 'outputType' is 'RegularSystem' or 'RS'        #
#       [r, cs] .......... 'outputType' is 'ConstructibleSet' or 'CS'     #
#   Where r is the rank of A for all parameter values that belonging to   #
#   the solution set of cs or rs. Together, all the constructible sets or #
#   regular systems form a partition of the input constructible set or    #
#   regular system.                                                       #
#                                                                         #
# EXAMPLE                                                                 #
#   > with(RegularChains):                                                #
#   > with(ParametricMatrixTools):                                        #
#   >                                                                     #
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
ComprehensiveRank := module()
    
    export ModuleApply;
    
    local
        init,
        init_F_H,
        init_rs,
        init_cs,
        
        processOptions,
        checkInput,
        
        implementation,
        
        convertToRS,
        
        # ALGORITHMS
        comprehensive_rank;
    
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
#    Same as ComprehensiveRank                                            #
# ----------------------------------------------------------------------- #
init := proc()

    local A, F, H, R, opts, cs, rs;

    # Check the number of arguments
    if nargs < 2 then
        error "Insufficient number of arguments";
    elif nargs > 5 then
        error "To many arguments";
    end if;

    if type(args[2], 'list') and type(args[3], 'list') then
        # ComprehensiveRank(A, F, H, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveRank called as ComprehensiveRank(A, F, H, R, options)");

        if nargs = 3 then
            error "Expected a fourth argument of a polynomial ring";
        end if;

        A := args[1];
        F := args[2];
        H := args[3];
        R := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_F_H(A, F, H, R, opts);

    elif type(args[2], 'list') then
        # ComprehensiveRank(A, F, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveRank called as ComprehensiveRank(A, F, R, options)");

        A := args[1];
        F := args[2];
        R := args[3];
        
        opts := processOptions({args[4..-1]});
        
        return init_F_H(A, F, [], R, opts);

    elif RC:-TRDis_regular_system(args[2]) then
        # ComprehensiveRank(A, rs, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveRank called as ComprehensiveRank(A, rs, R, options)");

        A  := args[1];
        rs := args[2];
        R  := args[3];
        
        opts := processOptions({args[4..-1]});
        
        return init_rs(A, rs, R, opts);

    elif RC:-TRDis_constructible_set(args[2]) then
        # ComprehensiveRank(A, cs, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveRank called as ComprehensiveRank(A, cs, R, options)");

        A  := args[1];
        cs := args[2];
        R  := args[3];
        
        opts := processOptions({args[4..-1]});
        
        return init_cs(A, cs, R, opts);
    
    elif RC:-TRDis_polynomial_ring(args[2]) then
        # ComprehensiveRank(A, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveRank called as ComprehensiveRank(A, R, options)");
        
        A  := args[1];
        R  := args[2];
        
        opts := processOptions({args[3..-1]});
        
        return init_F_H(A, [], [], R, opts);
        
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
#   See ComprehensiveRank header for specifications.                      #
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
#   F ...... List of polynomials representing equations                   #
#   H ...... List of polynomials representing inequations                 #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see ComprehensiveRank        #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveRank                                            #
# ----------------------------------------------------------------------- #
init_F_H := proc(A::Matrix, F::list(polynom), H::list(polynom), R::TRDring, opts::table, $)
    
    local p :: polynom,
          cs :: TRDcs;
    
    # Check the input for errors
    checkInput(A, R);
    
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
    
    return implementation(A, cs, R, opts);
    
end proc;


# ----------------------------------------------------------------------- #
# init_rs                                                                 #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   rs ..... Regular system                                               #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see ComprehensiveRank        #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveRank                                            #
# ----------------------------------------------------------------------- #
init_rs := proc(A::Matrix, rs::TRDrs, R::TRDring, opts::table, $)
    
    local cs :: TRDcs;
    
    # Check the input for errors
    checkInput(A, R);
    
    # rs must be a regular system
    if not RC:-TRDis_regular_system(rs, R) then
        error "Expected a regular system";
    end if;
    
    # Convert rs to a constructible set
    cs := RC_CST:-ConstructibleSet([rs], R);
    
    return implementation(A, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# init_cs                                                                 #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   cs ..... Constructible set                                            #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see ComprehensiveRank        #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveRank                                            #
# ----------------------------------------------------------------------- #
init_cs := proc(A::Matrix, cs::TRDcs, R::TRDring, opts, $)

    # Check the input for errors
    checkInput(A, R);

    # cs must be a constructible set
    if not RC:-TRDis_constructible_set(cs, R) then
        error "Expected a constructible set";
    end if;
    
    return implementation(A, cs, R, opts);

end proc;


# ----------------------------------------------------------------------- #
# checkInput                                                              #
#                                                                         #
# Checks for errors in the input values.                                  #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   R ...... Polynomial Ring                                              #
#   opts ... A table containing the options (see ComprehensiveRank        #
#            header)                                                      #
# ----------------------------------------------------------------------- #
checkInput := proc(A::Matrix, R::TRDring, $)
    
    local i :: posint, 
          j :: posint, 
          n :: nonnegint,
          m :: nonnegint;
    
    n, m := LA:-Dimension(A);
    
    # A must be at least a 1x1 matrix
    if n < 1 or m < 1 then
        error "Cannot compute rank of an empty matrix";
    end if;
    
    # All elements of A must be polynomials in R
    for i to n do
        for j to m do
            if not RC:-TRDis_poly(A[i,j], R) then
                error "Expected a matrix of polynomials in R";
            end if;
        end do;
    end do;
    
end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Computes the rank using the specified method and returns the specified  #
# type. Assume no errors in input values.                                 #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   cs ..... Constructible set                                            #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see ComprehensiveRank        #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveRank                                            #
# ----------------------------------------------------------------------- #
implementation := proc(AA::Matrix, cs::TRDcs, R::TRDring, opts::table, $)
    
    local rs :: TRDrs,
          rc :: TRDrc,
          A :: Matrix,
          r :: nonnegint,
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
            r := LA:-Rank(A);
            result := [op(result), [r, rs]];
            
        # Check if A is a zero matrix over rs
        elif isZeroMatrixOverRS(A, rs, R) then
            result := [op(result), [0, rs]];
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
        result := [op(result), op(comprehensive_rank(AA, csCompute, R))];
    end if;
    
    # Output options
    if opts['output_RS'] then
        result := convertToRS(result, R);
    end if;
    
    return result;
    
end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# EXTERNAL FILES
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveRank/comprehensive_rank.mpl>
$include <src/ComprehensiveRank/convertToRS.mpl>

end module;