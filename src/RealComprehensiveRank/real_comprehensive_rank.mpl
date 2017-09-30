# ======================================================================= #
# ======================================================================= #
#                                                                         #
# real_comprehensive_rank.mpl                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Sept. 29/2017                                               #
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
        constructTable,
        convertRankList,
        linearProjection_lrsas,
        linearProjection_rsas,
        SuggestVariableOrderWithParameters,
        add_list_polys_to_lrsas,
        getRealDimension;
    
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

    local n :: posint,
          m :: posint,
          y,
          lp,
          R_Union :: TRDring,
          lrsas_Int,
          rankTable,
          d,
          lrsas_singular,
          lrsas_nonsingular,
          out,
          i;
    
    # Get size of input matrix
    n, m := op(1,A);
    
    # Temporary variables for matrix
    y := seq('x'[i], i=1..m);
    
    # Get linear equations from matrix
    lp := convert(A.Vector([y]), list);
    
    # Use suggest variable order for improved performance
    R_Union := RC:-PolynomialRing(SuggestVariableOrderWithParameters(lp, R['variables']));
    
    # Pre-condition for case where A is singular
    d := LA:-Determinant(A);
    
    # Parameter values where A is singular
    lrsas_singular := add_list_polys_to_lrsas([d], lrsas, R);
    
    # If the matrix is non-singular for all parameter values
    # return the input conditions with full rank
    if nops(lrsas_singular) = 0 then
        return [[n, lrsas]];
    end if;
    
    # Parameter values where A is non-singular
    lrsas_nonsingular := RC_SAST:-Difference(lrsas, lrsas_singular, R);
    
    # Add the equations from the matrix
    lrsas_Int := add_list_polys_to_lrsas(lp, lrsas_singular, R_Union);
    
    # Call method to get rank
    rankTable := getRank(lrsas_Int, nops(R['variables']), R_Union);
    
    #return rankTable;
    out := convertRankList(rankTable);
    
    if nops(lrsas_nonsingular) <> 0 then
        out := [op(out), [n, lrsas_nonsingular]];
    end if;
    
    return out;
    
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

    local nVars, dim, j, i, rank, tab, lrsasNew, lrsas2;

    # Get the number of variables
    nVars := nops(R['variables']) - nParams;

    # Get dimension of each regular system
    dim := map(getRealDimension, lrsas, nVars, R);

    # arrange lrsas in order of least dimension to greatest dimension
    lrsasNew := NULL;
    for j from 0 to nVars do
        for i from 1 to nops(dim) do
            if evalb(dim[i] = j) then
                lrsasNew := lrsasNew, lrsas[i];
            end if;
        end do;
    end do;
    lrsasNew := [lrsasNew];

    # Compute Differences of all rsas with all rsas of greater of equal dimension
    lrsas2 := {};
    for i from 1 to nops(dim) - 1 do
        lrsas2 := lrsas2 union {op(RC_SAST:-Difference([lrsasNew[i]], lrsasNew[i+1..nops(lrsas)], R))};
    end do;
    lrsas2 := lrsas2 union {lrsasNew[-1]};
    lrsas2 := convert(lrsas2, list);

    # Construct the table with rank and regular_systems (pre-projection)
    tab := constructTable(lrsas2, nParams, R);

    # initialize rank table
    rank := table();
    for i from 0 to nVars do
        rank[i] := {};
    end do;

    # project variables onto parameter space
    for i from 0 to nVars do
        if nops(tab[i]) > 0 then
            rank[i] := [linearProjection_lrsas(convert(tab[i], list), nParams, R)];
        end if;
    end do;
    
    return rank;

end proc;


# ----------------------------------------------------------------------- #
# constructTable                                                          #
#                                                                         #
# Constructs a table such that table[r] corresponds to all regular        #
# semi-algebraic systems of rank r from the input list.                   #
#                                                                         #
# INPUT                                                                   #
#   lrsas ... List of regular_semi_algebraic_systems                      #
#   n ....... Number of parameters                                        #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#   Returns a table where table[r] contains a set of                      #
#   regular semi-algebraic systems corresponding to a rank of r.          #
# ----------------------------------------------------------------------- #
constructTable := proc(lrsas::TRDlrsas, n::nonnegint, R::TRDring, $)

    local nVars::nonnegint,
    tab :: table,
    rsas :: TRDrsas,
    d :: nonnegint,
    i :: nonnegint;

    # number of unknowns (x's)
    nVars := nops(R['variables']) - n;

    # initialize table to empty sets at each entry
    tab := table();
    for i from 0 to nVars do
        tab[i] := {};
    end do;

    # Compute the dimension of each rsas, add to table at
    # correct place
    for rsas in lrsas do
        d := getRealDimension(rsas, nVars, R);
        tab[nVars - d] := tab[nVars - d] union {rsas};
    end do;

    return tab;

end proc;


# ----------------------------------------------------------------------- #
# convertRankList                                                         #
#                                                                         #
# Converts a table where table[r] = lrsas where r is the computed rank    #
# and lrsas is a list of regular semi-algeraic systems to a list with     #
# elements of the form [r, lrsas].                                        #
#                                                                         #
# INPUT                                                                   #
#   rankTable ... Table to convert                                        #
#   R ........... Polynomial ring                                         #
#                                                                         #
# OUTPUT                                                                  #
#   List with elements indicating the rank                                #
# ----------------------------------------------------------------------- #
convertRankList := proc(rankTable::table, $)

    local idx :: list(nonnegint),
          i :: posint,
          outList :: list([nonnegint, TRDlrsas]);

    idx := sort([indices(rankTable,'nolist')]);

    outList := [];
    for i from 1 to nops(idx) do
        if nops(rankTable[idx[i]]) > 0 then
            outList := [op(outList), [idx[i], op(rankTable[idx[i]])]];
        end if;
    end do;
    
    return outList;
    
end proc;


# ----------------------------------------------------------------------- #
# getRealDimension                                                        #
#                                                                         #
# Computed the dimension (i.e. number of equations that do not contain    #
# the variables) in a regular semi-algebraic system.                      #
#                                                                         #
# INPUT                                                                   #
#   rsas ... Regular semi-algebraic system                                #
#   m ...... Number of unknowns                                           #
#   R ...... Polynomial ring                                              #
#                                                                         #
# OUTPUT                                                                  #
#   Returns a positive integer corresponding to the dimension of the      #
#   input rsas.                                                           #
# ----------------------------------------------------------------------- #
getRealDimension := proc(rsas, nVars, R, $)

    local x :: set(name),
          rc :: TRDrc,
          eqns :: list(polynom),
          eqn :: polynom,
          i :: nonnegint;
    
    # Get the variables
    x := {op((R['variables'])[1..nVars])};
    
    # Get the equations of rsas
    rc := RC_SAST:-RepresentingChain(rsas, R);
    eqns := RC:-Equations(rc, R);
    
    i := nVars;

    for eqn in eqns do
        if not evalb(nops(indets(eqn) intersect x) = 0) and not evalb(eqn = 0) then
            i := i - 1;
        end if;
    end do;

    return i;

end proc;


# ----------------------------------------------------------------------- #
# linearProjection_lrsas                                                  #
#                                                                         #
# Projects a list of regular semi-algebraic system onto the parameter     #
# space where all entries of lrsas are linear in the unknown variables.   #
#                                                                         #
# INPUT                                                                   #
#   rsas ... Regular semi-algebraic system                                #
#   n ...... Number of parameters                                         #
#   R ...... Polynomial ring                                              #
#                                                                         #
# OUTPUT                                                                  #
#   A lrsas without the equations containing the unknowns.                #
# ----------------------------------------------------------------------- #
linearProjection_lrsas := proc(lrsas, nParams, R)
    map(linearProjection_rsas, lrsas, nParams, R);
end proc;


# ----------------------------------------------------------------------- #
# linearProjection_rsas                                                   #
#                                                                         #
# Projects a regular semi-algebraic system onto the parameter space where #
# all entries of rsas are linear in the unknown variables.                #
#                                                                         #
# INPUT                                                                   #
#   rsas ... Regular semi-algebraic system                                #
#   n ...... Number of parameters                                         #
#   R ...... Polynomial ring                                              #
#                                                                         #
# OUTPUT                                                                  #
#   A rsas without the equations containing the unknowns.                 #
# ----------------------------------------------------------------------- #
linearProjection_rsas := proc(rsas::TRDrsas, nParams::nonnegint, R::TRDring, $)

    local R_params :: TRDring,
          nVars :: posint,
          ineqs :: list(polynom),
          qff :: TRDqff,
          rc :: TRDrc;
    
    nVars := nops(R['variables']) - nParams;
    
    # new polynomial ring without variables
    R_params := RC:-PolynomialRing(R['variables'][nVars+1..-1]);
    
    # list of inequalities
    ineqs := op(rsas)['PositiveInequalities'];
    
    # if there are no inequalities, make the empty list (Maple error!)
    if not type(ineqs, list) then
        ineqs := [];
    end if;

    qff := (op(rsas))['QuantifierFreeFormula'];
    rc := (op(rsas))['RegularChain'];
    
    rc := RC_CT:-Under(R['variables'][nVars], rc, R);

    return RC:-TRDmake_regular_semi_algebraic_system(qff, rc, ineqs, R_params);

end proc;


# ----------------------------------------------------------------------- #
# SuggestVariableOrderWithParameters                                      #
#                                                                         #
# Gives a variable ordering suggestion ensuring the parameters are less   #
# than any other variables.                                               #
#                                                                         #
# INPUT                                                                   #
#   lp ....... list of polynomials                                        #
#   params ... list of parameters                                         #
#                                                                         #
# OUTPUT                                                                  #
#   A list of variables.                                                  #
# ----------------------------------------------------------------------- #
SuggestVariableOrderWithParameters := proc(lp::list(polynom), params::list(name), $)
    
    local y :: list(name), v :: name;
    
    # Call standard suggest variable order
    y := RC:-SuggestVariableOrder(lp, params);
    
    # Make sure parameters are least variables
    for v in params do
        y := remove(x -> evalb(x=v), y);
    end do;
    
    return [op(LT:-Reverse(y)), op(params)];
    
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