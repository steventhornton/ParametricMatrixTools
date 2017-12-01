# ======================================================================= #
# ======================================================================= #
#                                                                         #
# real_comprehensive_rank.mpl                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 1/2017                                                 #
#                                                                         #
# Computes a complete case discussion of the rank of a matrix where the   #
# entries are multivariate polynomials whose indeterminants are treated   #
# as parameters. Computation is done modulo a regular system or           #
# constructible set. Parameters are assumed to be real valued and         #
# polynomial inequality constraints are allowed on the parameters.        #
#                                                                         #
# INPUT                                                                   #
#   A ....... Matrix                                                      #
#   lrsas ... Constructible set                                           #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [r_i, lrsas_i]                                                    #
#   where r_i is the rank of A for all parameter values that satisfy the  #
#   equations and inequations of all regular semi-algebraic systems in    #
#   lrsas_i. Together, all the lists of regular semi-algebraic systems    #
#   lrsas_i (for all i) represent a partition of the input list of        #
#   regular semi-algebraic systems.                                       #
#                                                                         #
# REFERENCES                                                              #
#                                                                         #
# ======================================================================= #
# ======================================================================= #
real_comprehensive_rank := module()

    export ModuleApply;

    local
        implementation_lrsas,
        getRank,
        convertRankTable,
        add_list_polys_to_lrsas;
    
    ModuleApply := proc(A::Matrix, lrsas::TRDlrsas, R::TRDring, $)
        return implementation_lrsas(A, lrsas, R);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ----------------------------------------------------------------------- #
# implementation_lrsas                                                    #
#                                                                         #
# Compute the rank of a parametric matrix.                                #
#                                                                         #
# INPUT/OUTPUT                                                            #
#    Same as real_comprehensive_rank                                      #
# ----------------------------------------------------------------------- #
implementation_lrsas := proc(A::Matrix, lrsas::TRDlrsas, R::TRDring, $)

    local m :: posint,
          y,
          i :: posint,
          lp :: list(polynom),
          yOrd :: list,
          R_Union :: TRDring,
          lrsas_linear :: TRDlrsas,
          rankTable :: table;
    
    # Get number of columns of A
    m := LA:-ColumnDimension(A);
    
    # Temporary variables for matrix
    y := seq('x'[i], i=1..m);
    
    # Get the equations from A
    lp := map(lhs, LA:-GenerateEquations(A, [y]));
    
    # Use SuggestVariableOrder huristics to pick the order of y
    yOrd := remove(a -> a in convert(R['variables'], set), RC:-SuggestVariableOrder(lp));
    
    # Join the parameters with the linear variables.
    R_Union := RC:-PolynomialRing([op(yOrd), op(R['variables'])]);
    
    # Add the equations from the matrix
    lrsas_linear := add_list_polys_to_lrsas(lp, lrsas, R_Union);
    
    # Compute the rank
    rankTable := getRank(lrsas_linear, nops(R['variables']), R_Union);
    
    return convertRankTable(rankTable, nops(yOrd));
    
end proc;


# ----------------------------------------------------------------------- #
# getRank                                                                 #
#                                                                         #
# Computes all possible rank values of a matrix with parametric           #
# polynomial entries with input algebraic polynomial constraints.         #
#                                                                         #
# INPUT                                                                   #
#   lrsas ..... List of regular semi-algebraic systems                    #
#   nParams ... Number of parameters                                      #
#   R ......... Polynomial ring (contains variables for both              #
#               matrix and parameters)                                    #
#                                                                         #
# OUTPUT                                                                  #
#   Returns a table where table[r] contains a set of regular              #
#   semi-algrbraic systems corresponding to a rank of r.                  #
# ----------------------------------------------------------------------- #
getRank := proc(lrsas, nParams, R::TRDring, $)

    local nVars :: posint,
          rankTable :: table,
          rsas :: TRDrsas,
          rc :: TRDrc,
          qff :: TRDqff,
          list_p :: list(polynom),
          dim :: nonnegint,
          i :: nonnegint,
          j :: posint;

    # Get the number of variables
    nVars := nops(R['variables']) - nParams;
    
    rankTable := table();
    
    for rsas in lrsas do
        
        # Extract the regular chain, positive inequalities, and
        rc := RC_SAST:-RepresentingChain(rsas, R);
        qff := RC_SAST:-RepresentingQuantifierFreeFormula(rsas);
        list_p := RC_SAST:-PositiveInequalities(rsas, R);
        
        # Compute the dimension of the regular chain
        dim := regularChainDimension(rc, nVars, R);
        
        # Remove any equations containing the linear variables
        rc := RC_CT:-Under(R['variables'][nVars], rc, R);
        
        # Reconstruct a regular semi-algebraic system
        rsas := RC:-TRDmake_regular_semi_algebraic_system(qff, rc, list_p, R);
        
        # Add to the rank table
        if not dim in [indices(rankTable, 'nolist')] then
            rankTable[dim] := [rsas];
        else
            rankTable[dim] := [op(rankTable[dim]), rsas];
        end if;
        
    end do;
    
    # Compute differences
    for i from 0 to nVars-1 do
        if not i in [indices(rankTable, 'nolist')] then next end if;
        for j from i+1 to nVars do
            if not j in [indices(rankTable, 'nolist')] then next end if;
            rankTable[i] := RC_SAST:-Difference(rankTable[i], rankTable[j], R);
        end do;
    end do;
    
    return(rankTable);
    
end proc;


# ----------------------------------------------------------------------- #
# convertRankTable                                                        #
#                                                                         #
# INPUT                                                                   #
#   t ... A table with non-negative integer indices and entries that are  #
#         lists or regular semi-algebraic systems.                        #
#   n ... Positive integer                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with entries                                                   #
#       [i, lrsas]                                                        #
#   where i is n-idx where idx is the index of t where lrsas occurs.      #
# ----------------------------------------------------------------------- #
convertRankTable := proc(t, n)

    local out,
    i;

    out := NULL;

    for i in [indices(t, 'nolist')] do
        if nops(t[i]) > 0 then
            out := out, [n - i, t[i]];
        end if;
    end do;

    return([out]);

end proc;


# ----------------------------------------------------------------------- #
# add_list_polys_to_lrsas                                                 #
#                                                                         #
# Adds the polynomials from a list (where each poly is = 0) to a list of  #
# regular semi-algebraic systems.                                         #
#                                                                         #
# INPUT                                                                   #
#   lp ...... List of polynomials                                         #
#   lrsas ... List of regular semi-algebraic systems                      #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#   A list of regular semi-algebraic systems :Intersection(V(lp), lrsas)  #
# ----------------------------------------------------------------------- #
add_list_polys_to_lrsas := proc(lp::{list(polynom),set(polynom)}, lrsas::TRDlrsas, R::TRDring, $)
    
    local lrsas_lp :: TRDlrsas;
    
    # What it needs to do, just faster!!
    lrsas_lp := RC:-RealTriangularize(lp,[],[],[],R);
    
    return RC_SAST:-Intersection(lrsas, lrsas_lp, R);
    
end proc;

end module;