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
real_comprehensive_rank := module()

    export ModuleApply;

    local
        implementation_lrsas,
        getRank,
        convertRankTable,
        add_list_polys_to_lrsas,
        getDimension;
    
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
        dim := getDimension(rc, nVars, R);
        
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



# 
# # ----------------------------------------------------------------------- #
# # constructTable                                                          #
# #                                                                         #
# # Constructs a table such that table[r] corresponds to all regular        #
# # semi-algebraic systems of rank r from the input list.                   #
# #                                                                         #
# # INPUT                                                                   #
# #   lrsas ... List of regular_semi_algebraic_systems                      #
# #   n ....... Number of parameters                                        #
# #   R ....... Polynomial ring                                             #
# #                                                                         #
# # OUTPUT                                                                  #
# #   Returns a table where table[r] contains a set of                      #
# #   regular semi-algebraic systems corresponding to a rank of r.          #
# # ----------------------------------------------------------------------- #
# constructTable := proc(lrsas::TRDlrsas, n::nonnegint, R::TRDring, $)
# 
#     local nVars::nonnegint,
#     tab :: table,
#     rsas :: TRDrsas,
#     d :: nonnegint,
#     i :: nonnegint;
# 
#     # number of unknowns (x's)
#     nVars := nops(R['variables']) - n;
# 
#     # initialize table to empty sets at each entry
#     tab := table();
#     for i from 0 to nVars do
#         tab[i] := {};
#     end do;
# 
#     # Compute the dimension of each rsas, add to table at
#     # correct place
#     for rsas in lrsas do
#         d := getRealDimension(rsas, nVars, R);
#         tab[nVars - d] := tab[nVars - d] union {rsas};
#     end do;
# 
#     return tab;
# 
# end proc;
# 
# 
# # ----------------------------------------------------------------------- #
# # convertRankList                                                         #
# #                                                                         #
# # Converts a table where table[r] = lrsas where r is the computed rank    #
# # and lrsas is a list of regular semi-algeraic systems to a list with     #
# # elements of the form [r, lrsas].                                        #
# #                                                                         #
# # INPUT                                                                   #
# #   rankTable ... Table to convert                                        #
# #   R ........... Polynomial ring                                         #
# #                                                                         #
# # OUTPUT                                                                  #
# #   List with elements indicating the rank                                #
# # ----------------------------------------------------------------------- #
# convertRankList := proc(rankTable::table, $)
# 
#     local idx :: list(nonnegint),
#           i :: posint,
#           outList :: list([nonnegint, TRDlrsas]);
# 
#     idx := sort([indices(rankTable,'nolist')]);
# 
#     outList := [];
#     for i from 1 to nops(idx) do
#         if nops(rankTable[idx[i]]) > 0 then
#             outList := [op(outList), [idx[i], op(rankTable[idx[i]])]];
#         end if;
#     end do;
#     
#     return outList;
#     
# end proc;
# 
# 
# # ----------------------------------------------------------------------- #
# # getRealDimension                                                        #
# #                                                                         #
# # Computed the dimension (i.e. number of equations that do not contain    #
# # the variables) in a regular semi-algebraic system.                      #
# #                                                                         #
# # INPUT                                                                   #
# #   rsas ... Regular semi-algebraic system                                #
# #   m ...... Number of unknowns                                           #
# #   R ...... Polynomial ring                                              #
# #                                                                         #
# # OUTPUT                                                                  #
# #   Returns a positive integer corresponding to the dimension of the      #
# #   input rsas.                                                           #
# # ----------------------------------------------------------------------- #
# # getRealDimension := proc(rsas, nVars, R, $)
# # 
# #     local x :: set(name),
# #           rc :: TRDrc,
# #           eqns :: list(polynom),
# #           eqn :: polynom,
# #           i :: nonnegint;
# #     
# #     # Get the variables
# #     x := {op((R['variables'])[1..nVars])};
# #     
# #     # Get the equations of rsas
# #     rc := RC_SAST:-RepresentingChain(rsas, R);
# #     eqns := RC:-Equations(rc, R);
# #     
# #     i := nVars;
# # 
# #     for eqn in eqns do
# #         if not evalb(nops(indets(eqn) intersect x) = 0) and not evalb(eqn = 0) then
# #             i := i - 1;
# #         end if;
# #     end do;
# # 
# #     return i;
# # 
# # end proc;
# 
# 
# # ----------------------------------------------------------------------- #
# # linearProjection_lrsas                                                  #
# #                                                                         #
# # Projects a list of regular semi-algebraic system onto the parameter     #
# # space where all entries of lrsas are linear in the unknown variables.   #
# #                                                                         #
# # INPUT                                                                   #
# #   rsas ... Regular semi-algebraic system                                #
# #   n ...... Number of parameters                                         #
# #   R ...... Polynomial ring                                              #
# #                                                                         #
# # OUTPUT                                                                  #
# #   A lrsas without the equations containing the unknowns.                #
# # ----------------------------------------------------------------------- #
# linearProjection_lrsas := proc(lrsas, nParams, R)
#     map(linearProjection_rsas, lrsas, nParams, R);
# end proc;
# 
# 
# # ----------------------------------------------------------------------- #
# # linearProjection_rsas                                                   #
# #                                                                         #
# # Projects a regular semi-algebraic system onto the parameter space where #
# # all entries of rsas are linear in the unknown variables.                #
# #                                                                         #
# # INPUT                                                                   #
# #   rsas ... Regular semi-algebraic system                                #
# #   n ...... Number of parameters                                         #
# #   R ...... Polynomial ring                                              #
# #                                                                         #
# # OUTPUT                                                                  #
# #   A rsas without the equations containing the unknowns.                 #
# # ----------------------------------------------------------------------- #
# linearProjection_rsas := proc(rsas::TRDrsas, nParams::nonnegint, R::TRDring, $)
# 
#     local R_params :: TRDring,
#           nVars :: posint,
#           ineqs :: list(polynom),
#           qff :: TRDqff,
#           rc :: TRDrc;
#     
#     nVars := nops(R['variables']) - nParams;
#     
#     # new polynomial ring without variables
#     R_params := RC:-PolynomialRing(R['variables'][nVars+1..-1]);
#     
#     # list of inequalities
#     ineqs := op(rsas)['PositiveInequalities'];
#     
#     # if there are no inequalities, make the empty list (Maple error!)
#     if not type(ineqs, list) then
#         ineqs := [];
#     end if;
# 
#     qff := (op(rsas))['QuantifierFreeFormula'];
#     rc := (op(rsas))['RegularChain'];
#     
#     rc := RC_CT:-Under(R['variables'][nVars], rc, R);
# 
#     return RC:-TRDmake_regular_semi_algebraic_system(qff, rc, ineqs, R_params);
# 
# end proc;


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


# ----------------------------------------------------------------------- #
# getDimension                                                            #
#                                                                         #
# Compute the dimension of a regular chain. The dimension is              #
# nVars - (# of equations containing one or more of the largest nVars     #
# variables of R). The input regular chain is assumed to be linear in     #
# each of the nVars largest variables of R.                               #
#                                                                         #
# INPUT                                                                   #
#   rc ...... Regular chain                                               #
#   nVars ... Number of variables                                         #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#   Returns a positive integer corresponding to the dimension of the      #
#   input regular chain.                                                  #
# ----------------------------------------------------------------------- #
getDimension := proc(rc::TRDrc, nVars::nonnegint, R::TRDring, $)

    local x :: set(name),
          eqns :: list(polynom),
          eqn :: polynom,
          i :: nonnegint;

    # Get the largest 'nVars' variables or R
    x := {op((R['variables'])[1..nVars])};

    # Get the equations
    eqns := RC:-Equations(rc, R);

    # Initialize the dimension at nVars
    i := nVars;

    for eqn in eqns do
        if not evalb(nops(indets(eqn) intersect x) = 0) and not evalb(eqn = 0) then
            i := i - 1;
        end if;
    end do;

    return i;

end proc;

end module;