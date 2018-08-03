# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ctd_comprehensive_rank.mpl                                              #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Aug. 3/2018                                                 #
#                                                                         #
# Computes a complete case discussion of the rank of a matrix where the   #
# entries are multivariate polynomials whose indeterminants are treated   #
# as parameters. Computation is done modulo a regular system or           #
# constructible set.                                                      #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [r_i, cs_i]                                                       #
#   where r_i is the rank of A for all parameter values that satisfy the  #
#   equations and inequations of cs_i. Together, all the constructible    #
#   sets cs_i (for all values of i) represent partition of the input      #
#   constructible set.                                                    #
#                                                                         #
# REFERENCES                                                              #
#                                                                         #
# ======================================================================= #
# ======================================================================= #
ctd_comprehensive_rank := proc(A::Matrix, cs::TRDcs, R::TRDring, $)
    
    local n, d, X, x, linear_eqs, R2, lrc, lcs, lrc_rank, X_vec, lrs, yOrd, cs_S, i;
    
    # Get the size of A
    n := LA:-ColumnDimension(A);
    
    # Get the number of parameters in the polynomial ring
    d := nops(R['variables']);
    
    # Generate a vector of variables for the linear system
    X := [seq(x[i], i = 1..n)];
    X_vec := Vector(X);
    
    # Convert matrix to linear equations in X
    linear_eqs := convert(A.X_vec, list);
    
    # Use SuggestVariableOrder huristics to pick the order of y
    yOrd := remove(a -> a in convert(R['variables'], set), RC:-SuggestVariableOrder(linear_eqs));
    yOrd := LT:-Reverse(yOrd);
    
    # Join the parameters with the linear variables.
    R2 := RC:-PolynomialRing([op(yOrd), op(R['variables'])]);
    
    # Compute a constructible set of the linear system
    cs_S := RC_CST:-GeneralConstruct(linear_eqs, [], R2);
    
    # Intersect the input constraints on the parameters with the linear system
    cs_S := RC_CST:-Intersection(cs_S, cs, R2);
    
    # Compute the comprehensive triangular decomposition
    lrs, lcs := RC_PST:-ComprehensiveTriangularize(cs_S, d, R2);
    
    # Extract the regular chains from each regular system in lrs
    lrc := map(rs -> RC_CST:-RepresentingChain(rs, R), lrs);
    
    # Compute the "rank" of each regular chain in lrc (i.e. count the number
    # of equations containing any of the variables in X)
    lrc_rank := map(rc -> nops(RC_CT:-Upper(R['variables'][1], rc, R2)), lrc);
    
    # Get the minimum rank of the regular systems (regular chains) corresponding
    # to each constructible set
    return map(a -> [min(lrc_rank[a[2]]), a[1]], lcs);
    
end proc;
