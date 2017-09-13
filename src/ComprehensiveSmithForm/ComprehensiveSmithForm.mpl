# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ComprehensiveSmithForm.mpl                                              #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Sept. 12/2017                                               #
#                                                                         #
# Computes the Smith normal form of a matrix where the entries are        #
# parametric univariate polynomials. Computation is done modulo a         #
# regular system or constructible set.                                    #
#                                                                         #
# CALLING SEQUENCE                                                        #
#    ComprehensiveSmithForm(A, v, R, options)                             #
#    ComprehensiveSmithForm(A, v, rs, R, options)                         #
#    ComprehensiveSmithForm(A, v, cs, R, options)                         #
#    ComprehensiveSmithForm(A, v, F, R, options)                          #
#    ComprehensiveSmithForm(A, v, F, H, R, options)                       #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   v .... Variable                                                       #
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
#   outputMatrices ... Default = S                                        #
#                      S, U, V, [S, U, V], [U, V], [S, U], [S, V]         #
#                         - Which matrices to return (see below)          #
#                                                                         #
#   lazy ............. Default = false                                    #
#                      If true, only one branch of the computation is     #
#                      returned (significatly faster). Otherwise, all     #
#                      branches are returned.                             #
#                                                                         #
#   algorithm ........ Default = minors                                   #
#                      Specify which algorithm to use (Only the minors    #
#                      algorithm is currently implemented).               #
#                                                                         #
# OPTION COMPATIBILITY                                                    #
#   - Only the Smith form can be returned when the minors alogrithm is    #
#     used, the transformation matrices are not computed with this        #
#     algorithm.                                                          #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements in one of the following forms:                   #
#       [S, rs] .......... 'outputType' is 'RegularSystem' or 'RS' and    #
#                                  'outputMatrices' = 'S'                 #
#       [S, cs] ......... 'outputType' is 'ConstructibleSet' or 'CS' and  #
#                                  'outputMatrices' = 'S'                 #
#       [U, rs] ......... 'outputType' is 'RegularSystem' or 'RS' and     #
#                                  'outputMatrices' = 'U'                 #
#       [U, cs] ......... 'outputType' is 'ConstructibleSet' or 'CS' and  #
#                                  'outputMatrices' = 'U'                 #
#       [V, rs] ......... 'outputType' is 'RegularSystem' or 'RS' and     #
#                                  'outputMatrices' = 'V'                 #
#       [V, cs] ......... 'outputType' is 'ConstructibleSet' or 'CS' and  #
#                                  'outputMatrices' = 'V'                 #
#       [S, U, rs] ...... 'outputType' is 'RegularSystem' or 'RS' and     #
#                                  'outputMatrices' = {'S', 'U'}          #
#       [S, U, cs] ...... 'outputType' is 'ConstructibleSet' or 'CS' and  #
#                                  'outputMatrices' = {'S', 'U'}          #
#       [S, V, rs] ...... 'outputType' is 'RegularSystem' or 'RS' and     #
#                                  'outputMatrices' = {'S', 'U', 'V'}     #
#       [S, V, cs] ...... 'outputType' is 'ConstructibleSet' or 'CS' and  #
#                                  'outputMatrices' = {'S', 'V'}          #
#       [S, U, V, rs] ... 'outputType' is 'RegularSystem' or 'RS' and     #
#                                  'outputMatrices' = {'S', 'U', 'V'}     #
#       [S, U, V, cs] ... 'outputType' is 'ConstructibleSet' or 'CS' and  #
#                                  'outputMatrices' = {'S', 'U', 'V'}     #
#       [U, V, rs] ...... 'outputType' is 'RegularSystem' or 'RS' and     #
#                                  'outputMatrices' = {'U', 'V'}          #
#       [U, V, cs] ...... 'outputType' is 'ConstructibleSet' or 'CS' and  #
#                                  'outputMatrices' = {'U', 'V'}          #
#   Where S is the Smith Normal Form of A for all parameter values that   #
#   satisfy the equations and inequations of cs or rs. U and V are the    #
#   left and right transformation matrices respectively such that         #
#       S = U.A.V                                                         #
#   Together, all the constructible sets or regular systems form a        #
#   partition of the input constructible set or regular system.           #
#                                                                         #
# ASSUMPTIONS                                                             #
#   - v is the greatest variable in R.                                    #
#   - If lists or polynomials (F and H) are input, v can not be a         #
#     variable in any polynomial in these lists.                          #
#   - If a constructible set or regular system is input, v cannot be a    #
#     variable in any equations or inequations.                           #
#                                                                         #
# EXAMPLE                                                                 #
#   > with(RegularChains):                                                #
#   > with(ParametricMatrixTools):                                        #
#   >                                                                     #
#   > R := PolynomialRing([x,a]):                                         #
#   > A := Matrix([[ 0,  a,  1,  1,  1],                                  #
#                  [ 2, -2,  0, -2, -4],                                  #
#                  [ 0,  0,  1,  1,  3],                                  #
#                  [-6,  0, -3, -1, -3],                                  #
#                  [ 2,  2,  2,  2,  4]]):                                #
#   > Ax := x*IdentityMatrix(5) - A:                                      #
#   >                                                                     #
#   > result := ComprehensiveSmithForm(Ax, x, [], R):                     #
#   >                                                                     #
#   > map(factor, result[1][1]), Display(result[1][2], R);                #
#         [1,0,0,0,0]                                                     #
#         [0,1,0,0,0]                              {a-4 <> 0              #
#         [0,0,1,0,0]                              {a-1 <> 0              #
#         [0,0,0,1,0]                              {23a^2-61a+68 <> 0     #
#         [0,0,0,0,-(x-2)(x^2+2)(-x^2+2a-4)]                              #
#   >                                                                     #
#   > map(factor, result[2][1]), Display(result[2][2], R);                #
#         [1,0,0,0,0]                                                     #
#         [0,1,0,0,0]                                                     #
#         [0,0,1,0,0]                              {23a^2-61a+68 = 0      #
#         [0,0,0,1,0]                                                     #
#         [0,0,0,0,-(x-2)(x^2+2)(-x^2+2a-4)]                              #
#   >                                                                     #
#   > map(factor, result[3][1]), Display(result[3][2], R);                #
#         [1,0,0,0,  0]                                                   #
#         [0,1,0,0,  0]                                                   #
#         [0,0,1,0,  0]                     {a-4 = 0                      #
#         [0,0,0,x-2,0]                                                   #
#         [0,0,0,0,  (x-2)*(x+2)*(x^2+2)]                                 #
#   >                                                                     #
#   > map(factor, result[4][1]), Display(result[4][2], R);                #
#         [1,0,0,0,    0]                                                 #
#         [0,1,0,0,    0]                                                 #
#         [0,0,1,0,    0]                     {a-1 = 0                    #
#         [0,0,0,x^2+2,0]                                                 #
#         [0,0,0,0,    (x-2)*(x^2+2)]                                     #
#                                                                         #
# REFERENCES                                                              #
#                                                                         #
# TO DO                                                                   #
#   1. Correct ordering of output matrices based on input option          #
#   2. Implement outputType option                                        #
#   3. Implement lazy option                                              #
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
ComprehensiveSmithForm := module()
    
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
        comprehensive_smith_form_minors;
    
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
#    Same as ComprehensiveSmithForm                                       #
# ----------------------------------------------------------------------- #
init := proc()

    local A, v, F, H, R, opts, cs, rs;

    # Check the number of arguments
    if nargs < 3 then
        error "Insufficient number of arguments";
    elif nargs > 9 then
        error "To many arguments";
    end if;

    if type(args[3], 'list') and type(args[4], 'list') then
        # ComprehensiveSmithForm(A, v, F, H, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveSmithForm called as ComprehensiveSmithForm(A, v, F, H, R, options)");

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
        # ComprehensiveSmithForm(A, v, F, R, options)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveSmithForm called as ComprehensiveSmithForm(A, v, F, R, options)");

        A := args[1];
        v := args[2];
        F := args[3];
        R := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_F_H(A, v, F, [], R, opts);

    elif RC:-TRDis_regular_system(args[3]) then
        # ComprehensiveSmithForm(A, v, rs, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveSmithForm called as ComprehensiveSmithForm(A, v, rs, R, options)");

        A  := args[1];
        v  := args[2];
        rs := args[3];
        R  := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_rs(A, v, rs, R, opts);

    elif RC:-TRDis_constructible_set(args[3]) then
        # ComprehensiveSmithForm(A, v, cs, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveSmithForm called as ComprehensiveSmithForm(A, v, cs, R, options)");

        A  := args[1];
        v  := args[2];
        cs := args[3];
        R  := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_cs(A, v, cs, R, opts);
    
    elif RC:-TRDis_polynomial_ring(args[3]) then
        # ComprehensiveSmithForm(A, v, R, opts)
        userinfo(2, 'ParametricMatrixTools', "ComprehensiveSmithForm called as ComprehensiveSmithForm(A, v, R, options)");
        
        A  := args[1];
        v  := args[2];
        R  := args[3];
        
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
#       output_S                                                          #
#       output_U                                                          #
#       output_V                                                          #
#       output_S_ord                                                      #
#       output_U_ord                                                      #
#       output_V_ord                                                      #
#       lazy                                                              #
#       algorithm_minors                                                  #
#       algorithm_elementary                                              #
#   See ComprehensiveSmithForm header for specifications.                 #
# ----------------------------------------------------------------------- #
processOptions := proc(opts_in::set(equation), $) :: table;

    local opts::table(),
          tab_opts,
          opt::equation,
          opt_name,
          opt_value,
          i;
    
    # Default values
    opts['output_CS'] := true;
    opts['output_RS'] := false;
    opts['output_S']  := true;
    opts['output_U']  := false;
    opts['output_V']  := false;
    opts['output_S_ord'] := 1;
    opts['output_U_ord'] := 0;
    opts['output_V_ord'] := 0;
    opts['lazy']      := false;
    opts['algorithm_minors'] := true;
    opts['algorithm_elementary'] := false;
    
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
                    if opt_value[i] = 'S' then
                        opts['output_S'] := true;
                        opts['output_S_ord'] := i;
                    elif opt_value[i] = 'U' then
                        opts['output_U'] := true;
                        opts['output_U_ord'] := i;
                    elif opt_value[i] = 'V' then
                        opts['output_V'] := true;
                        opts['output_V_ord'] := i;
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
                if opt_value = 'minors' then
                    opts['algorithm_minors'] := true;
                    opts['algorithm_elementary'] := false;
                elif opt_value = 'elementary' then
                    opts['algorithm_minors'] := false;
                    opts['algorithm_elementary'] := true;
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
    if opts['algorithm_minors'] then
        if opts['output_U'] or opts['output_V'] then
            error "Ouput matrices U and V are incompatible with the minors algorithm.";
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
#   v ...... Variable                                                     #
#   F ...... List of polynomials representing equations                   #
#   H ...... List of polynomials representing inequations                 #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see ComprehensiveSmithForm   #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveSmithForm                                       #
# ----------------------------------------------------------------------- #
init_F_H := proc(A::Matrix(square), v::name, F::list(polynom), H::list(polynom), R::TRDring, opts::table, $)
    
    local p :: polynom,
          cs :: TRDcs;
    
    # Check the input for errors
    checkInput(A, v, R);
    
    # All elements of F must be polynomials in R
    for p in F do
        if not RC:-TRDis_poly(p, R) then
            error "Invalid polynomial in F";
        end if;
        if v in indets(p) then
            error "Invalid polynomial in F";
        end if;
    end do;
    
    # All elements of H must be polynomials in R
    for p in H do
        if not RC:-TRDis_poly(p, R) then
            error "Invalid polynomial in H";
        end if;
        if v in indets(p) then
            error "Invalid polynomial in F";
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
#   opts ... A table containing the options (see ComprehensiveSmithForm   #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveSmithForm                                       #
# ----------------------------------------------------------------------- #
init_rs := proc(A::Matrix(square), v::name, rs::TRDrs, R::TRDring, opts::table, $)
    
    local cs :: TRDcs;
    
    # Check the input for errors
    checkInput(A, v, R);
    
    # rs must be a regular system
    if not RC:-TRDis_regular_system(rs, R) then
        error "Expected a regular system";
    end if;
    
    if not isUnder(rs, v, R) then
        error "Input regular system must not contain conditions on %1", v;
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
#   opts ... A table containing the options (see ComprehensiveSmithForm   #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveSmithForm                                       #
# ----------------------------------------------------------------------- #
init_cs := proc(A::Matrix(square), v::name, cs::TRDcs, R::TRDring, opts, $)

    # Check the input for errors
    checkInput(A, v, R);

    # cs must be a constructible set
    if not RC:-TRDis_constructible_set(cs, R) then
        error "Expected a constructible set";
    end if;
    
    if not isUnder(cs, v, R) then
        error "Input constructible set must not contain conditions on %1", v;
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
#   R ...... Polynomial Ring                                              #
#   opts ... A table containing the options (see ComprehensiveSmithForm   #
#            header)                                                      #
# ----------------------------------------------------------------------- #
checkInput := proc(A::Matrix(square), v::name, R::TRDring, $)
    
    local i :: posint, 
          j :: posint, 
          n :: nonnegint,
          m :: nonnegint;
    
    # v must be the greatest variable of R
    if not isGreatestVariable(v, R) then
        error "v must be the greatest variable of R";
    end if;
    
    n, m := LA:-Dimension(A);
    
    # A must be at least a 1x1 matrix
    if n < 1 then
        error "Cannont compute Smith Form of an empty matrix";
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
# Computes the Smith form using the specified method and returns the      #
# specified type. Assume no errors in input values.                       #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   v ...... Variable                                                     #
#   cs ..... Constructible set                                            #
#   R ...... Polynomial ring                                              #
#   opts ... A table containing the options (see ComprehensiveSmithForm   #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ComprehensiveSmithForm                                       #
# ----------------------------------------------------------------------- #
implementation := proc(AA::Matrix(square), v::name, cs::TRDcs, R::TRDring, opts::table, $)
    
    local rs :: TRDrs,
          rc :: TRDrc,
          A :: 'Matrix'(square),
          S :: 'Matrix'(square),
          U :: 'Matrix'(square),
          V :: 'Matrix'(square),
          result := [],
          lrsCompute :: TRDlrs := [],
          csCompute :: TRDcs;
          
    
    # Check for zero or constant matrices for each regular system in cs
    for rs in RC_CST:-RepresentingRegularSystems(cs, R) do
        
        # Map SparsePseudoRemainder to A
        rc := RC_CST:-RepresentingChain(rs, R);
        A := map(RC:-SparsePseudoRemainder, AA, rc, R);
        
        # Check if A is either a constant matrix, or has v as the only 
        # indeterminant
        if isConstantMatrix(A, v, R) then
            
            if opts['output_S'] then
                S := LA:-SmithForm(A, v, 'output'='S');
                result := [op(result), [S, rs]];
            elif opts['output_U'] then
                U := LA:-SmithForm(A, v, 'output'='U');
                result := [op(result), [U, rs]];
            elif opts['output_V'] then
                V := LA:-SmithForm(A, v, 'output'='V');
                result := [op(result), [V, rs]];
            elif opts['output_S'] and opts['output_U'] then
                S, U := LA:-SmithForm(A, v, 'output'=['S','U']);
                result := [op(result), [S, U, rs]];
            elif opts['output_S'] and opts['output_V'] then
                S, V := LA:-SmithForm(A, v, 'output'=['S','V']);
                result := [op(result), [S, V, rs]];
            elif opts['output_S'] and opts['output_U'] and opts['output_V'] then
                S, U, V := LA:-SmithForm(A, v, 'output'=['S','U','V']);
                result := [op(result), [S, U, V, rs]];
            elif opts['output_U'] and opts['output_V'] then
                U, V := LA:-SmithForm(A, v, 'output'=['U', 'V']);
                result := [op(result), [U, V, rs]];
            end if;
        
        # Check if A is a zero matrix over rs
        elif isZeroMatrixOverRS(A, rs, R) then
            
            if opts['output_S'] then
                result := [op(result), [A, rs]];
            elif opts['output_U'] then
                result := [op(result), [LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['output_V'] then
                result := [op(result), [LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['output_S'] and opts['output_U'] then
                result := [op(result), [A, LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['output_S'] and opts['output_V'] then
                result := [op(result), [A, LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['output_S'] and opts['output_U'] and opts['output_V'] then
                result := [op(result), [A, LA:-IdentityMatrix(LA:-Dimension(A)), LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['output_U'] and opts['output_V'] then
                result := [op(result), [LA:-IdentityMatrix(LA:-Dimension(A)), LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
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
    if opts['algorithm_elementary'] then
        error "The ComprehensiveSmithForm elementry algorithm has not been implemented yet.";
    elif opts['algorithm_minors'] then
        result := [op(result), op(comprehensive_smith_form_minors(AA, v, csCompute, R))];
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
$include <src/ComprehensiveSmithForm/comprehensive_smith_form_minors.mpl>
$include <src/ComprehensiveSmithForm/convertToRS.mpl>

end module;