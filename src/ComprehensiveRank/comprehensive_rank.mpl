# ======================================================================= #
# ======================================================================= #
#                                                                         #
# comprehensive_rank.mpl                                                  #
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
comprehensive_rank := module()

    export ModuleApply;

    local
        implementation_cs,
        getRank,
        convertRankTable,
        add_list_polys_to_cs;
    
    ModuleApply := proc(A::Matrix, cs::TRDcs, R::TRDring, $)
        return implementation_cs(A, cs, R);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ----------------------------------------------------------------------- #
# implementation_cs                                                       #
#                                                                         #
# Compute the rank of a parametric matrix.                                #
#                                                                         #
# INPUT/OUTPUT                                                            #
#    Same as comprehensive_rank                                           #
# ----------------------------------------------------------------------- #
implementation_cs := proc(A::Matrix, cs::TRDcs, R::TRDring, $)

    local m :: posint,
          y,
          i :: posint,
          lp :: list(polynom),
          yOrd :: list,
          R_Union :: TRDring,
          cs_linear :: TRDcs,
          lrs :: TRDlrs,
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
    
    # Add the equations from the matrix to the input constructible set
    cs_linear := add_list_polys_to_cs(lp, cs, R_Union);
    
    # Convert to a list of regular systems
    lrs := RC:-TRDregular_systems(cs_linear, R_Union);
    
    # Compute the rank
    rankTable := getRank(lrs, nops(R['variables']), R_Union);
    
    return convertRankTable(rankTable, nops(yOrd));
    
end proc;


# ----------------------------------------------------------------------- #
# convertRankTable                                                        #
#                                                                         #
# INPUT                                                                   #
#   t ... A table with non-negative integer indices and entries that are  #
#         constructible sets.                                             #
#   n ... Positive integer                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with entries                                                   #
#       [i, cs]                                                           #
#   where i is n-idx where idx is the index of t where cs occurs.         #
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
# getRank                                                                 #
#                                                                         #
# Counts the number of of variables in the input regular system...?       #
#                                                                         #
# INPUT                                                                   #
#   lrs ....... List of regular systems                                   #
#   nParams ... number of parameters                                      #
#   R ......... Polynomial Ring (contains variables for both matrix       #
#               and parameters)                                           #
#                                                                         #
# OUTPUT                                                                  #
#   Returns a table where table[r] contains a set of regular_systems      #
#   corresponding to a rank of r.                                         #
# ----------------------------------------------------------------------- #
getRank := proc(lrs::TRDlrs, nParams::posint, R::TRDring, $)

    local nVars :: posint,
          rankTable :: table,
          rs :: TRDrs,
          H :: list(polynom),
          T :: TRDrc,
          dim :: nonnegint,
          i :: nonnegint,
          j :: nonnegint;
          
    # Get the number of variables
    nVars := nops(R['variables']) - nParams;
    
    rankTable := table();
    
    for rs in lrs do
        # Extract the regular chain and set of inequations from the regular 
        # system
        H := RC_CST:-RepresentingInequations(rs, R);
        T := RC_CST:-RepresentingChain(rs, R);
        
        # Compute the dimension of the regular chain
        dim := regularChainDimension(T, nVars, R);
        
        # Remove any equations containing the linear variables
        T := RC_CT:-Under(R['variables'][nVars], T, R);
        
        # Reconstruct a regular system
        rs := RC_CST:-RegularSystem(T, H, R);
        
        # Add to the rank table
        if not dim in [indices(rankTable, 'nolist')] then
            rankTable[dim] := {rs}
        else
            rankTable[dim] := `union`(rankTable[dim], {rs});
        end if;
        
    end do;
    
    # Convert each set in rankTable to a constructible set
    for i in [indices(rankTable, 'nolist')] do
        rankTable[i] := RC_CST:-ConstructibleSet(convert(rankTable[i], 'list'), R);
    end do;
    
    # Compute differences
    for i from 0 to nVars-1 do
        if not i in [indices(rankTable, 'nolist')] then next end if;
        for j from i+1 to nVars do
            if not j in [indices(rankTable, 'nolist')] then next end if;
            rankTable[i] := RC:-TRDAlgebraicDifference(rankTable[i], rankTable[j], R);
        end do;
    end do;
    
    return(rankTable);
    
end proc;


# ----------------------------------------------------------------------- #
# add_list_polys_to_cs                                                    #
#                                                                         #
# Adds the polynomials from a list (where each poly is = 0)               #
# to a constructible set.                                                 #
#                                                                         #
# INPUT                                                                   #
#   lp ... List of polynomials                                            #
#   CS ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   Intersection(V(lp), CS)                                               #
# ----------------------------------------------------------------------- #
add_list_polys_to_cs := proc(lp::{list,set}, CS, R)

    local lrs, rs, lcs, rc, ieqs;

    lrs := RC_CST:-RepresentingRegularSystems(CS, R);
    
    lcs := NULL;
    for rs in lrs do
        rc := RC_CST:-RepresentingChain(rs, R);
        ieqs := RC_CST:-RepresentingInequations(rs, R);
        lcs := lcs, RC_CST:-GeneralConstruct(lp, rc, ieqs, R);
    end do;

    return ListUnion([lcs], R);

end proc;

end module;