# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ComprehensiveJordanForm.mpl                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 10/2017                                                #
#                                                                         #
# Computes the Jordan canonical form of a matrix where the entries are    #
# mutlivariate polynomials. Computation is done modulo a regular system   #
# or constructible set.                                                   #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   ComprehensiveJordanForm(A, R, options)                                #
#   ComprehensiveJordanForm(A, rs, R, options)                            #
#   ComprehensiveJordanForm(A, cs, R, options)                            #
#   ComprehensiveJordanForm(A, F, R, options)                             #
#   ComprehensiveJordanForm(A, F, H, R, options)                          #
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
#   outputMatrices ... Default = J                                        #
#                      J, Q, [J, Q]                                       #
#                         - Which matrices to return (see below)          #
#                                                                         #
#   lazy ............. Default = true                                     #
#                      If true, only one branch of the computation is     #
#                      returned (significatly faster). Otherwise, all     #
#                      branches are returned.                             #
#                                                                         #
#   algorithm ........ Default = snf_minors                               #
#                      Specify which algorithm to use (Only the           #
#                      snf_minors algorithm is currently implemented).    #
#                                                                         #
# OPTION COMPATIBILITY                                                    #
#   - Only the Jordan form can be returned when the snf_minors alogrithm  #
#     is used, the transformation matrix is not computed with this        #
#     algorithm.                                                          #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements in one of the following forms:                   #
#       [J, rs] .......... 'outputType' is 'RegularSystem' or 'RS' and    #
#                                  'outputMatrices' = 'J'                 #
#       [J, cs] ......... 'outputType' is 'ConstructibleSet' or 'CS' and  #
#                                  'outputMatrices' = 'J'                 #
#       [Q, rs] ......... 'outputType' is 'RegularSystem' or 'RS' and     #
#                                  'outputMatrices' = 'Q'                 #
#       [Q, cs] ......... 'outputType' is 'ConstructibleSet' or 'CS' and  #
#                                  'outputMatrices' = 'Q'                 #
#       [J, Q, rs] ...... 'outputType' is 'RegularSystem' or 'RS' and     #
#                                  'outputMatrices' = ['J', 'Q']          #
#       [J, Q, cs] ...... 'outputType' is 'ConstructibleSet' or 'CS' and  #
#                                  'outputMatrices' = ['J', 'Q']          #
#   Where J is the Jordan canonical form of A for all parameter values    #
#   that satisfy the equations and inequations of cs or rs. Q is the      #
#   similarity transformation matrix such that                            #
#       J = Q^(-1) A Q                                                    #
#   Together, all the constructible sets or regular systems form a        #
#   partition of the input constructible set or regular system.           #
#                                                                         #
# ASSUMPTIONS                                                             #
#                                                                         #
# EXAMPLE                                                                 #
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
ComprehensiveJordanForm := module()

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
        comprehensive_jordan_form_snf_minors;

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
#    Same as ComprehensiveJordanForm                                      #
# ----------------------------------------------------------------------- #
init := proc()

    local A :: 'Matrix'(square),
          F :: list(polynom),
          H :: list(polynom),
          R :: TRDring,
          rs :: TRDrs,
          cs :: TRDcs,
          opts :: table;

    if nargs < 2 then
        error "Insufficient number of arguments";
    elif nargs > 8 then
        error "To many arguments";
    end if;

    if type(args[2], 'list') and type(args[3], 'list') then
        # ComprehensiveJordanForm(A, F, H, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveJordanForm called as ComprehensiveJordanForm(A, F, H, R, options)");

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
        # ComprehensiveJordanForm(A, F, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveJordanForm called as ComprehensiveJordanForm(A, F, R, options)");

        A := args[1];
        F := args[2];
        R := args[3];

        opts := processOptions({args[4..-1]});

        return init_F_H(A, F, [], R, opts);

    elif RC:-TRDis_regular_system(args[2]) then
        # ComprehensiveJordanForm(A, rs, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveJordanForm called as ComprehensiveJordanForm(A, rs, R, options)");

        A  := args[1];
        rs := args[2];
        R  := args[3];

        opts := processOptions({args[4..-1]});

        return init_rs(A, rs, R, opts);

    elif RC:-TRDis_constructible_set(args[2]) then
        # ComprehensiveJordanForm(A, cs, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveJordanForm called as ComprehensiveJordanForm(A, cs, R, options)");

        A  := args[1];
        cs := args[2];
        R  := args[3];

        opts := processOptions({args[4..-1]});

        return init_cs(A, cs, R, opts);
    
    elif RC:-TRDis_polynomial_ring(args[2]) then
        # ComprehensiveJordanForm(A, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveJordanForm called as ComprehensiveJordanForm(A, R, options)");
        
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
#       output_F                                                          #
#       output_Q                                                          #
#       output_F_ord                                                      #
#       output_Q_ord                                                      #
#       lazy                                                              #
#       algorithm_snf_minors                                              #
#       algorithm_snf_elementary                                          #
#   See ComprehensiveJordanForm header for specifications.                #
# ----------------------------------------------------------------------- #
processOptions := proc(opts_in::set(equation), $) :: table;

    local opts::table,
          tab_opts,
          opt::equation,
          opt_name,
          opt_value,
          i;

    # Default values
    opts['output_CS'] := true;
    opts['output_RS'] := false;
    opts['output_J']  := true;
    opts['output_Q']  := false;
    opts['output_J_ord'] := 1;
    opts['output_Q_ord'] := 0;
    opts['lazy']      := false;
    opts['algorithm_snf_minors'] := true;
    opts['algorithm_snf_elementary'] := false;

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
            elif opt_name = ('outputMatrices') then
                opt_value := rhs(opt);
                opt_value := `if`(type(opt_value, 'list'), opt_value, [opt_value]);
                for i to nops(opt_value) do
                    if opt_value[i] = 'F' then
                        opts['output_J'] := true;
                        opts['output_J_ord'] := i;
                    elif opt_value[i] = 'Q' then
                        opts['output_Q'] := true;
                        opts['output_Q_ord'] := i;
                    else
                        error "incorrect option value: %1", opt;
                    end if;
                end do;
            elif opt_name = ('lazy') then
                opt_value := rhs(opt);
                if not type(opts['lazy'], 'truefalse') then
                    error "incorrect option value: %1", opt;
                end if;
                opts['lazy'] := opt_value;
            elif opt_name = ('algorithm') then
                opt_value := rhs(opt);
                if opt_value = 'snf_minors' then
                    opts['algorithm_snf_minors'] := true;
                    opts['algorithm_snf_elementary'] := false;
                elif opt_value = 'snf_elementary' then
                    opts['algorithm_snf_minors'] := false;
                    opts['algorithm_snf_elementary'] := true;
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

    # Check option compatibility
    if opts['algorithm_snf_minors'] or opts['algorithm_snf_elementary'] then
        if opts['output_Q'] then
            error "Ouput matrix Q is incompatible with the snf_minors and snf_elementary algorithms.";
        end if;
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
#   A ...... Matrix of multivariate polynomials                           #
#   F ...... List of polynomials representing equations                   #
#   H ...... List of polynomials representing inequations                 #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see                          #
#            ComprehensiveJordanForm header)                              #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveJordanForm                                      #
# ----------------------------------------------------------------------- #
init_F_H := proc(A::Matrix(square), F::list(polynom), H::list(polynom), R::TRDring, opts::table, $)

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
#   opts ... A table containing the options (see                          #
#            ComprehensiveJordanForm header)                              #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveJordanForm                                      #
# ----------------------------------------------------------------------- #
init_rs := proc(A::Matrix(square), rs::TRDrs, R::TRDring, opts::table, $)

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
#   opts ... A table containing the options (see                          #
#            ComprehensiveJordanForm header)                              #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveJordanForm                                      #
# ----------------------------------------------------------------------- #
init_cs := proc(A::Matrix(square), cs::TRDcs, R::TRDring, opts, $)

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
# ----------------------------------------------------------------------- #
checkInput := proc(A::Matrix(square), R::TRDring, $)

    local i :: posint, 
          j :: posint, 
          n :: nonnegint,
          m :: nonnegint;

    n, m := LA:-Dimension(A);

    # A must be at least a 1x1 matrix
    if n < 1 then
        error "Cannont compute Jordan form of an empty matrix";
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
# Computes the Jordan form using the specified method and returns the     #
# specified type. Assume no errors in input values.                       #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   cs ..... Constructible set                                            #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see                          #
#            ComprehensiveJordanForm header)                              #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveJordanForm                                      #
# ----------------------------------------------------------------------- #
implementation := proc(AA::Matrix(square), cs::TRDcs, R::TRDring, opts::table, $)

    local rs :: TRDrs,
          rc :: TRDrc,
          A :: 'Matrix'(square),
          J :: 'Matrix'(square),
          Q :: 'Matrix'(square),
          result := [],
          lrsCompute :: TRDlrs := [],
          csCompute :: TRDcs;


    # Check for zero or constant matrices for each regular system in cs
    for rs in RC_CST:-RepresentingRegularSystems(cs, R) do
        
        # Map SparsePseudoRemainder to A
        rc := RC_CST:-RepresentingChain(rs, R);
        A := map(RC:-SparsePseudoRemainder, AA, rc, R);

        # Check if A is either a constant matrix
        if isConstantMatrix(A, R) then

            if opts['output_J'] then
                Q := ParametricMatrixTools:-JordanForm(A, 'output'='J');
                result := [op(result), [J, rs]];
            elif opts['output_Q'] then
                Q := ParametricMatrixTools:-JordanForm(A, 'output'='J');
                result := [op(result), [Q, rs]];
            elif opts['output_J'] and opts['output_Q'] then
                J, Q := ParametricMatrixTools:-JordanForm(A, 'output'=['J','Q']);
                result := [op(result), [J, Q, rs]];
            end if;

        # Check if A is a zero matrix over rs
        elif isZeroMatrixOverRS(A, rs, R) then

            if opts['output_J'] then
                result := [op(result), [Matrix(LA:-Dimension(A)), rs]];
            elif opts['output_Q'] then
                result := [op(result), [LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['output_J'] and opts['output_Q'] then
                result := [op(result), [A, LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            end if;
        else
            lrsCompute := [op(lrsCompute), rs];
        end if;    
        
    end do;

    # Convert all regular systems in result list to constructible sets
    result := map(input -> [RC_CST:-ConstructibleSet([input[-1]], R), 
                            input[1..-2]], result);
    
    # Convert lrsCompute to a constructible set
    csCompute := RC_CST:-ConstructibleSet(lrsCompute, R);
    
    # Call the correct algorithm
    if opts['algorithm_snf_elementary'] then
        error "The ComprehensiveJordanForm snf_elementry algorithm has not been implemented yet.";
    elif opts['algorithm_snf_minors'] then
        result := [op(result), op(comprehensive_jordan_form_snf_minors(AA, csCompute, R))];
    end if;

    return result;

end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# EXTERNAL FILES
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveJordanForm/comprehensive_jordan_form_snf_minors.mpl>

end module;