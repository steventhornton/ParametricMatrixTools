# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ParametricSmithForm_Minors.mpl                                          #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 24/2016                                                #
#                                                                         #
# Computes the Smith Normal Form of a matrix where the entries are        #
# univariate polynomials with multivariate coefficients in the            #
# parameters. Computation is done modulo a regular system or              #
# constructible set.                                                      #
#                                                                         #
# CALLING SEQUENCE                                                        #
#    ParametricSmithForm(A, v, rs, R, options)                            #
#    ParametricSmithForm(A, v, cs, R, options)                            #
#    ParametricSmithForm(A, v, F, R, options)                             #
#    ParametricSmithForm(A, v, F, H, R, options)                          #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   v .... Variable                                                       #
#   rs ... Regular system                                                 #
#   cs ... Constructible set                                              #
#   F .... List of polynomials over R representing equations              #
#   H .... List of polynomials over R representing inequations            #
#   R .... Polynomial Ring                                                #
#                                                                         #
# OPTIONS                                                                 #
#   outputType ....... Default = CS                                       #
#                      ConstructibleSet or CS:                           #
#                         - Output will contain constructible sets        #
#                      RegularSystem or RS:                              #
#                           - Output will contain regular systems         #
#                                                                         #
#   outputMatrices ... Default = S                                        #
#                      S, U, V, [S, U, V], [U, V], [S, U], [S, V]         #
#                         - Which matrices to return (see below)          #
#                                                                         #
#   lazy ............. Default = true                                     #
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
ParametricSmithForm := module()
    
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
        ParametricSmithForm_Minors;
        #ParametricSmithForm_Elementary,
        
        # OTHER METHODS
        #clean_rs,
        #clean_cs,
        #convertToRS,
        #convertToCS;
    
    ModuleApply := proc() :: {list([Matrix, TRDrs]), list([Matrix, TRDcs]), list([Matrix, Matrix, TRDrs]), list([Matrix, Matrix, TRDcs]), list([Matrix, Matrix, Matrix, TRDrs]), list([Matrix, Matrix, Matrix, TRDcs])};
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
#    Same as ParametricSmithForm                                          #
# ----------------------------------------------------------------------- #
init := proc()
    
    local A :: 'Matrix'(square),
          v :: name,
          F :: list(polynom),
          H :: list(polynom),
          R :: TRDring,
          rs :: TRDrs,
          cs :: TRDcs,
          opts :: table;
    
    if nargs < 4 then
        error "Insufficient number of arguments";
    elif nargs > 9 then
        error "To many arguments";
    end if;
    
    if type(args[3], 'list') and type(args[4], 'list') then
        # ParametricSmithForm(A, v, F, H, R, opts)
        
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
        # ParametricSmithForm(A, v, F, R, opts)
        
        A := args[1];
        v := args[2];
        F := args[3];
        R := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_F_H(A, v, F, [], R, opts);
        
    elif RC:-TRDis_regular_system(args[3]) then
        # ParametricSmithForm(A, v, rs, R, opts)
        
        A  := args[1];
        v  := args[2];
        rs := args[3];
        R  := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_rs(A, v, rs, R, opts);
        
    elif RC:-TRDis_constructible_set(args[3]) then
        # ParametricSmithForm(A, v, cs, R, opts)
        
        A  := args[1];
        v  := args[2];
        cs := args[3];
        R  := args[4];
        
        opts := processOptions({args[5..-1]});
        
        return init_cs(A, v, cs, R, opts);
        
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
#   A table with indices                                                  #
#       'outputType'                                                      #
#       'outputMatrices'                                                  #
#       'lazy'                                                            #
#       'algorithm'                                                       #
#   See ParametricSmithForm header for specifications.                    #
# ----------------------------------------------------------------------- #
processOptions := proc(opts_in::set(equation), $) :: table;

    local opts::table(),
          opt::equation;

    # Default values
    opts['outputType']     := 'CS';
    opts['outputMatrices'] := 'S';
    opts['lazy']           := true;
    opts['algorithm']      := 'Elementary';

    # Process each option
    for opt in opts_in do
        if lhs(opt) in {indices(opts, 'nolist')} then
            opts[lhs(opt)] := rhs(opt);
        else
            error "%1 is not a valid option", lhs(opt);
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
#   opts ... A table containing the options (see ParametricSmithForm      #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ParametricSmithForm                                          #
# ----------------------------------------------------------------------- #
init_F_H := proc(A::Matrix(square), v::name, F::list(polynom), H::list(polynom), R::TRDring, opts::table, $)
    
    local p :: polynom,
          cs :: TRDcs;
    
    # Check the input for errors
    checkInput(A, v, R, opts);
    
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
#   opts ... A table containing the options (see ParametricSmithForm      #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ParametricSmithForm                                          #
# ----------------------------------------------------------------------- #
init_rs := proc(A::Matrix(square), v::name, rs::TRDrs, R::TRDring, opts::table, $)
    
    local cs :: TRDcs;
    
    # Check the input for errors
    checkInput(A, v, R, opts);
    
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
#   opts ... A table containing the options (see ParametricSmithForm      #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ParametricSmithForm                                          #
# ----------------------------------------------------------------------- #
init_cs := proc(A::Matrix(square), v::name, cs::TRDcs, R::TRDring, opts, $)

    # Check the input for errors
    checkInput(A, v, R, opts);

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
#   opts ... A table containing the options (see ParametricSmithForm      #
#            header)                                                      #
# ----------------------------------------------------------------------- #
checkInput := proc(A::Matrix(square), v::name, R::TRDring, opts::table, $)
    
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
    
    # Check the outputMatrices option
    if not opts['outputMatrices'] in {'S', 'U', 'V', {'U', 'V'}, {'S', 'U'}, {'S', 'V'}, {'S', 'U', 'V'}} then
        error "outputMatrices option must be either S, U, V, {U, V}, {S, U}, {S, V}, {S, U, V}";
    end if;
    
    # Check the outputType option
    if not opts['outputType'] in {'RegularSystem', 'RS', 'ConstructibleSet', 'CS'} then
        error "output must be either RegularSystem, RS, ConstructibleSet or CS";
    end if;
    
    # Check the algorithm option
    if not opts['algorithm'] in {'Minors', 'Elementary'} then
        error "Algorithm option must be either Minors or Elementary";
    end if;
    
    # Check the lazy option
    if not type(opts['lazy'], 'truefalse') then
        error "lazy option must be a boolean values";
    end if;
    
    # Check for option incompatability
    if (opts['outputMatrices'] in {'U', 'V', {'U', 'V'}, {'S', 'U'}, {'S', 'V'}, {'S', 'U', 'V'}}) and opts['algorithm'] in {'Minors'} then
        error "Output options algorithm=Minors and U or V in outputMatrices are not compatible";
    end if;
    
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
#   opts ... A table containing the options (see ParametricSmithForm      #
#            header)                                                      #
#                                                                         #
# OUTPUT                                                                  #
#    Same as ParametricSmithForm                                          #
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
        if isConstantMatrix(A, v) then
            
            if opts['outputMatrices'] = "S" then
                S := LA:-SmithForm(A, v, 'output'='S');
                result := [op(result), [S, rs]];
            elif opts['outputMatrices'] = "U" then
                U := LA:-SmithForm(A, v, 'output'='U');
                result := [op(result), [U, rs]];
            elif opts['outputMatrices'] = "V" then
                V := LA:-SmithForm(A, v, 'output'='V');
                result := [op(result), [V, rs]];
            elif opts['outputMatrices'] = {"S", "U"} then
                S, U := LA:-SmithForm(A, v, 'output'=['S','U']);
                result := [op(result), [S, U, rs]];
            elif opts['outputMatrices'] = {"S", "V"} then
                S, V := LA:-SmithForm(A, v, 'output'=['S','V']);
                result := [op(result), [S, V, rs]];
            elif opts['outputMatrices'] = {"S", "U", "V"} then
                S, U, V := LA:-SmithForm(A, v, 'output'=['S','U','V']);
                result := [op(result), [S, U, V, rs]];
            elif opts['outputMatrices'] = {"U", "V"} then
                U, V := LA:-SmithForm(A, v, 'output'=['U', 'V']);
                result := [op(result), [U, V, rs]];
            end if;
        
        # Check if A is a zero matrix over rs
        elif isZeroMatrix(A, rs, R) then
            
            if opts['outputMatrices'] = "S" then
                result := [op(result), [A, rs]];
            elif opts['outputMatrices'] = "U" then
                result := [op(result), [LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['outputMatrices'] = "V" then
                result := [op(result), [LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['outputMatrices'] = {"S", "U"} then
                result := [op(result), [A, LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['outputMatrices'] = {"S", "V"} then
                result := [op(result), [A, LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['outputMatrices'] = {"S", "U", "V"} then
                result := [op(result), [A, LA:-IdentityMatrix(LA:-Dimension(A)), LA:-IdentityMatrix(LA:-Dimension(A)), rs]];
            elif opts['outputMatrices'] = {"U", "V"} then
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
    if opts['algorithm'] = 'Elementary' then
        #SList := ParametricSmithForm_Elementary(AA, v, csCompute, R, opts);
        printf("The ParametricSmithForm Elementry algorithm has not been implemented yet.\n");
    elif opts['algorithm'] = 'Minors' then
        result := [op(result), op(ParametricSmithForm_Minors(AA, v, csCompute, R))];
    end if;
    
    # DELETE THIS
    return result;
    
    # Convert result list into the correct type for output
    # if opts['outputType'] in {'ConstructibleSet', 'CS'} then
    #     result := clean_cs(result, R);
    # else
    #     result := convertToRS(result, R);
    #     result := clean_rs(result, R);
    # end if;
    # 
    # return result;
    
end proc:


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# OTHER METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ------------------------------------------------------------------------------
# clean_cs
#
# Combine any elements of the output list that are the same.
#
# INPUT
#    l ... A list with elements of the form
#            [A1, cs], [A1, A2, cs] or [A1, A2, A3, cs]
#          where A1, A2 and A3 are matrices and cs is a constructible set
#    R ... Polynomial Ring
#
# OUTPUT
#    A list of the same form as the input, might have fewer elements if
#     combination happened.
# ------------------------------------------------------------------------------
# clean_cs := proc(l, R::table, $)
#     
#     local ll::{list([Matrix, table]), list([Matrix, Matrix, table]), list([Matrix, Matrix, Matrix, table])},
#           item::{[Matrix, table], [Matrix, Matrix, table], [Matrix, Matrix, Matrix, table]},
#           csTasks::Stack := SimpleStack(),
#           rs_i::table,
#           A_i::{[Matrix], [Matrix, Matrix], [Matrix, Matrix, Matrix]},
#           out_lrs := [],
#           equalAll::truefalse,
#           noMatch::truefalse,
#           A_j::{[Matrix], [Matrix, Matrix], [Matrix, Matrix, Matrix]},
#           lrs_j::list(table),
#           rs_j::table,
#           j::posint,
#           k::posint,
#           cs::table,
#           out::{list([Matrix, table]), list([Matrix, Matrix, table]), list([Matrix, Matrix, Matrix, table])} := [];
#           
#     
#     # Convert all constructible sets to regular systems
#     ll := convertToRS(l, R);
#     
#     # Add everything from ll to the stack
#     for item in ll do
#         csTasks:-push(item)
#     end do;
#     
#     while not csTasks:-empty() do
#         
#         item := csTasks:-pop();
#         rs_i := item[-1];
#         A_i := item[1..-2];
#         
#         noMatch := true;
#         
#         for j to nops(out_lrs) while noMatch do
#             
#             lrs_j := out_lrs[j][-1];
#             A_j   := out_lrs[j][1..-2];
#             
#             equalAll := true;
#             
#             for rs_j in lrs_j do
#                 for k to nops(A_i) while equalAll do
#                     if not isZeroMatrixOverRS(A_i[k] - A_j[k], rs_i, R) or
#                        not isZeroMatrixOverRS(A_i[k] - A_j[k], rs_j, R) then
#                         equalAll := false;
#                     end if;
#                 end do;
#                 
#                 if not equalAll then
#                     break;
#                 end if;
#             end do;
#             
#             if equalAll then
#                 noMatch := false;
#                 out_lrs[j] := [op(A_j), [op(lrs_j), rs_i]];
#             end if;
#             
#         end do;
#         
#         if noMatch then
#             out_lrs := [op(out_lrs), [op(A_i), [rs_i]]];
#         end if;
#         
#     end do;
#     
#     # Convert lists of regular systems in out list to constructible sets
#     for j to nops(out_lrs) do
#         
#         lrs_j := out_lrs[j][-1];
#         A_j := out_lrs[j][1..-2];
#         
#         cs := RC_CST:-ConstructibleSet(lrs_j, R);
#     
#         out := [op(out), [op(A_j), cs]];
#     
#     end do;
#     
#     return out;
#     
# end proc:


# ------------------------------------------------------------------------------
# clean_rs
#
# Call SparsePseudoRemainder on each matrix
#
# INPUT
#    l ... A list with elements of the form
#            [A1, rs], [A1, A2, rs] or [A1, A2, A3, rs]
#          where A1, A2 and A3 are matrices and rs is a regular system
#    R ... Polynomial Ring
#
# OUTPUT
#    A list of the same form as the input, might have fewer elements if
#     combination happened.
# ------------------------------------------------------------------------------
# clean_rs := proc(l::{list([Matrix, table]), list([Matrix, Matrix, table]), list([Matrix, Matrix, Matrix, table])}, R::table, $) :: {list([Matrix, table]), list([Matrix, Matrix, table]), list([Matrix, Matrix, Matrix, table])};
# 
#     local item::{[Matrix, table], [Matrix, Matrix, table]},
#           A::{[Matrix], [Matrix, Matrix], [Matrix, Matrix, Matrix]},
#           i::posint,
#           rs::table,
#           rc::table,
#           out::{list([Matrix, table]), list([Matrix, Matrix, table])}:=[];
# 
#     for item in l do
#         
#         rs := item[-1];
#         A  := item[1..-2];
#         
#         rc := RC_CST:-RepresentingChain(rs, R);
#         
#         for i to nops(A) do
#             A[i] := map(RC:-SparsePseudoRemainder, A[i], rc, R)
#         end do;
#         
#         out := [op(out), [op(A), rs]];
#         
#     end do;
#     
#     return out;
# 
# end proc:


# ------------------------------------------------------------------------------
# convertToRS
#
# Split all constructible sets up into lists or regular systems
#
# INPUT
#    l ... A list with elements of the form
#            [A1, cs], [A1, A2, cs] or [A1, A2, A3, cs]
#          where A1, A2 and A3 are matrices and cs is a constructible set
#    R ... Polynomial Ring
#
# OUTPUT
#    A list with elements of the form
#        [A1, cs], [A1, A2, cs] or [A1, A2, A3, cs]
#    where A1, A2 and A3 are matrices and rs is a regular system.
# ------------------------------------------------------------------------------
# convertToRS := proc(l:: {list([Matrix, table]), list([Matrix, Matrix, table]), list([Matrix, Matrix, Matrix, table])}, R::table, $) :: {list([Matrix, table]), list([Matrix, Matrix, table]), list([Matrix, Matrix, Matrix, table])};
# 
#     local out::{list([Matrix, table]), list([Matrix, Matrix, table]), list([Matrix, Matrix, Matrix, table])} := [],
#           rs::table,
#           item::{[Matrix, table], [Matrix, Matrix, table], [Matrix, Matrix, Matrix, table]},
#           A::{[Matrix], [Matrix, Matrix], [Matrix, Matrix, Matrix]},
#           cs::table;
#     
#     for item in l do
#         
#         cs := item[-1];
#         A := item[1..-2];
#         
#         for rs in RC_CST:-RepresentingRegularSystems(cs, R) do
#             out := [op(out), [op(A), rs]];
#         end do;
#         
#     end do;
#     
#     return out;
# 
# end proc:


# ------------------------------------------------------------------------------
# convertToCS
#
# Convert all regular systems in list to constructible sets.
#
# INPUT
#    l ... A list with elements of the form
#            [A1, rs], [A1, A2, rs] or [A1, A2, A3, rs]
#          where A1, A2 and A3 are matrices and rs is a regular system
#    R ... Polynomial Ring
#
# OUTPUT
#    A list with elements of the form
#        [A1, cs], [A1, A2, cs] or [A1, A2, A3, cs]
#    where A1, A2 and A3 are matrices and cs is a constructible system consisting 
#    of one regular system. The output list has the same number of elements as 
#    the input list.
# ------------------------------------------------------------------------------
# convertToCS := proc(l:: {list([Matrix, table]), list([Matrix, Matrix, table]), list([Matrix, Matrix, Matrix, table])}, R::table, $) :: {list([Matrix, table]), list([Matrix, Matrix, table]), list([Matrix, Matrix, Matrix, table])};
#     
#     local out::{list([Matrix, table]), list([Matrix, Matrix, table]), list([Matrix, Matrix, Matrix, table])} := [],
#           rs::table,
#           item::{[Matrix, table], [Matrix, Matrix, table], [Matrix, Matrix, Matrix, table]},
#           A::{[Matrix], [Matrix, Matrix], [Matrix, Matrix, Matrix]},
#           cs::table;
#     
#     # Determine if elements are of the form [F, rs] or [F, Q, rs]
#     for item in l do
#         rs := item[-1];        # The regular system
#         A := item[1..-2];    # List with the matrices
#         
#         cs := RC_CST:-ConstructibleSet([rs], R);
#         
#         out := [op(out), [op(A), cs]];
#         
#     end do;
# 
#     return out;
# 
# end proc:


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# EXTERNAL FILES
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <SmithForm/ParametricSmithForm_Minors.mpl>

end module: