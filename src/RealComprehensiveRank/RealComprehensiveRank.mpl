# ======================================================================= #
# ======================================================================= #
#                                                                         #
# RealComprehensiveRank.mpl                                               #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 24/2017                                                #
#                                                                         #
# Computes a complete case discussion of the rank of a matrix where the   #
# entries are multivariate polynomials whose indeterminants are treated   #
# as parameters. Parameters are assumed to be real valued and computation #
# is done over polynomial equality, inequation and inequality (strict     #
# and non-strict) constraints on the parameters.                          #
#                                                                         #
# CALLING SEQUENCE                                                        #
#    RealComprehensiveRank(A, R)                                          #
#    RealComprehensiveRank(A, rsas, R)                                    #
#    RealComprehensiveRank(A, lrsas, R)                                   #
#    RealComprehensiveRank(A, F, R)                                       #
#    RealComprehensiveRank(A, F, H, R)                                    #
#    RealComprehensiveRank(A, F, N, P, H, R)                              #
#                                                                         #
# INPUT                                                                   #
#   A ....... Matrix                                                      #
#   rsas .... Regular semi-algebraic system                               #
#   lrsas ... List of regular semi-algebraic systems                      #
#   F ....... List of polynomials over R representing equations           #
#   N ....... List or set of polynomials (non-negativity)                 #
#   P ....... List or set of polynomials (positivity)                     #
#   H ....... List of polynomials over R representing inequations         #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements in one of the form:                              #
#       [r, lrsas]                                                        #
#   Where r is the rank of A for all parameter values that belonging to   #
#   the solution set of lrsas. Together, all the lists of regular         #
#   semi-algebraic systems form a partition of the input constraints on   #
#   the parameters.                                                       #
#                                                                         #
# EXAMPLE                                                                 #
#                                                                         #
# REFERENCES                                                              #
#                                                                         #
# ======================================================================= #
# ======================================================================= #
RealComprehensiveRank := module()
    
    export ModuleApply;
    
    local
        init,
        init_F_N_P_H,
        init_lrsas,
        
        checkInput,
        
        implementation,
        
        # ALGORITHMS
        real_comprehensive_rank;
    
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
#    Same as RealComprehensiveRank                                        #
# ----------------------------------------------------------------------- #
init := proc()

    local A, F, N, P, H, R, lrsas;
    
    # Check the number of arguments
    if nargs < 2 then
        error "Insufficient number of arguments";
    elif nargs > 6 then
        error "To many arguments";
    end if;
    
    if nargs = 2 then
        # RealComprehensiveRank(A, R)
        userinfo(2, 'ParametricMatrixTools', "RealComprehensiveRank called as RealComprehensiveRank(A, R)");
        
        A := args[1];
        F := [];
        N := [];
        P := [];
        H := [];
        R := args[2];
        
        return init_F_N_P_H(A, F, N, P, H, R);
        
    elif nargs = 3 and type(args[2], 'TRDlrsas')  then
        # RealComprehensiveRank(A, lrsas, R)
        userinfo(2, 'ParametricMatrixTools', "RealComprehensiveRank called as RealComprehensiveRank(A, lrsas, R)");
        
        A := args[1];
        lrsas := args[2];
        R := args[3];
        
        return init_lrsas(A, lrsas, R);
        
    elif nargs = 3 and type(args[2], 'TRDrsas') then
        # RealComprehensiveRank(A, rsas, R)
        userinfo(2, 'ParametricMatrixTools', "RealComprehensiveRank called as RealComprehensiveRank(A, rsas, R)");
        
        A := args[1];
        lrsas := [args[2]];
        R := args[3];
        
        return init_lrsas(A, lrsas, R);
        
    elif nargs = 3 and  type(args[2], list(polynom)) then
        # RealComprehensiveRank(A, F, R)
        userinfo(2, 'ParametricMatrixTools', "RealComprehensiveRank called as RealComprehensiveRank(A, F, R)");
        
        A := args[1];
        F := args[2];
        N := [];
        P := [];
        H := [];
        R := args[3];
        
        return init_F_N_P_H(A, F, N, P, H, R);
        
    elif nargs = 4 then
        # RealComprehensiveRank(A, F, H, R)
        userinfo(2, 'ParametricMatrixTools', "RealComprehensiveRank called as RealComprehensiveRank(A, F, H, R)");
        
        A := args[1];
        F := args[2];
        N := [];
        P := [];
        H := args[3];
        R := args[4];
        
        return init_F_N_P_H(A, F, N, P, H, R);
        
    elif nargs = 6 then
        # RealComprehensiveRank(A, F, N, P, H, R)
        userinfo(2, 'ParametricMatrixTools', "RealComprehensiveRank called as RealComprehensiveRank(A, F, N, P, H, R)");
        
        A := args[1];
        F := args[2];
        N := args[3];
        P := args[4];
        H := args[5];
        R := args[6];
        
        return init_F_N_P_H(A, F, N, P, H, R);
        
    else
        error "Invalid arguments";
    end if;

end proc;


# ----------------------------------------------------------------------- #
# init_F_N_P_H                                                            #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   A ...... Matrix of multivariate polynomials                           #
#   F ...... List of polynomials over R representing equations            #
#   N ...... List or set of polynomials (non-negativity)                  #
#   P ...... List or set of polynomials (positivity)                      #
#   H ...... List of polynomials over R representing inequations          #
#   R ...... Polynomial ring                                              #
#                                                                         #
# OUTPUT                                                                  #
#    Same as RealComprehensiveRank                                        #
# ----------------------------------------------------------------------- #
init_F_N_P_H := proc(A::Matrix, F::list(polynom), N::list(polynom), P::list(polynom), H::list(polynom), R::TRDring, $)
    
    local p :: polynom,
          lrsas :: TRDlrsas;
    
    # Check the input for errors
    checkInput(A, R);
    
    # All elements of F must be polynomials in R
    for p in F do
        if not RC:-TRDis_poly(p, R) then
            error "Invalid polynomial in F";
        end if;
    end do;
    
    # All elements of N must be polynomials in R
    for p in N do
        if not RC:-TRDis_poly(p, R) then
            error "Invalid polynomial in N";
        end if;
    end do;
    
    # All elements of P must be polynomials in R
    for p in P do
        if not RC:-TRDis_poly(p, R) then
            error "Invalid polynomial in P";
        end if;
    end do;
    
    # All elements of H must be polynomials in R
    for p in H do
        if not RC:-TRDis_poly(p, R) then
            error "Invalid polynomial in H";
        end if;
    end do;
    
    # Convert F, N, P, and H to a list of regular semi-algebraic systems
    lrsas := RC:-RealTriangularize(F, N, P, H, R, 'output'='list');
    
    return implementation(A, lrsas, R);
    
end proc;


# ----------------------------------------------------------------------- #
# init_lrsas                                                              #
#                                                                         #
# Checks for errors in the input values, calls appropriate implementation #
# if all checks pass.                                                     #
#                                                                         #
# INPUT                                                                   #
#   A ........ Matrix of multivariate polynomials                         #
#   lrsas ... List of regular semi-algebraic systems                      #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#    Same as RealComprehensiveRank                                        #
# ----------------------------------------------------------------------- #
init_lrsas := proc(A::Matrix, lrsas::TRDlrsas, R::TRDring, $)
    
    local rsas :: TRDrsas;
    
    # Check the input for errors
    checkInput(A, R);
    
    # lrsas must be a list of regular semi-algebraic systems
    for rsas in lrsas do
        if not RC:-TRDis_regular_semi_algebraic_system(rsas, R) then
            error "Expected a list of regular semi-algebraic system";
        end if;
    end do;
    
    return implementation(A, lrsas, R);

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
#   A ....... Matrix of multivariate polynomials                          #
#   lrsas ... List of regular semi-algebraic systems                      #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#    Same as RealComprehensiveRank                                        #
# ----------------------------------------------------------------------- #
implementation := proc(AA::Matrix, lrsas::TRDlrsas, R::TRDring, $)
    
    local rsas :: TRDrsas,
          rc :: TRDrc,
          A :: Matrix,
          r :: nonnegint,
          result := [],
          lrsasCompute :: TRDlrsas := [];
    
    # Check for zero or constant matrices for each regular system in cs
    for rsas in lrsas do
        
        # Map SparsePseudoRemainder to A
        rc := RC_SAST:-RepresentingChain(rsas, R);
        A := map(RC:-SparsePseudoRemainder, AA, rc, R);
        
        # Check if A is a constant matrix
        if isConstantMatrix(A, R) then
            r := LA:-Rank(A);
            result := [op(result), [r, rsas]];
            
        # Check if A is a zero matrix over rsas
        elif isZeroMatrixOverRSAS(A, rsas, R) then
            result := [op(result), [0, rsas]];
        else
            lrsasCompute := [op(lrsasCompute), rsas];
        end if;    
        
    end do;
    
    # Call the algorithm
    if nops(lrsasCompute) > 0 then
        result := [op(result), op(real_comprehensive_rank(AA, lrsasCompute, R))];
    end if;
    
    return result;
    
end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# EXTERNAL FILES
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/RealComprehensiveRank/real_comprehensive_rank.mpl>

end module;