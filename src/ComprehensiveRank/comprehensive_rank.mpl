# ======================================================================= #
# ======================================================================= #
#                                                                         #
# comprehensive_rank.mpl                                                  #
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
        constructTable,
        convertRankList,
        linearProjection,
        add_list_polys_to_cs,
        getDimension;
    
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
          lp :: list(polynom),
          R_Union :: TRDring,
          cs_linear :: TRDcs,
          lrs :: TRDlrs,
          rankTable :: table,
          out :: list([nonnegint, TRDcs]),
          i :: posint;

    # Get number of columns of A
    m := LA:-ColumnDimension(A);

    # Temporary variables for matrix
    y := seq('x'[i], i=1..m);
    
    # Get linear equations from matrix
    lp := convert(A.Vector([y]), list);

    # Join the parameters with the new variables.
    R_Union := RC:-PolynomialRing([y, op(R['variables'])]);
    
    # Add the equations from the matrix
    cs_linear := add_list_polys_to_cs(lp, cs, R_Union);

    # Convert to list of regular systems
    lrs := RC:-TRDregular_systems(cs_linear, R_Union);

    # Call method to get rank
    rankTable := getRank(lrs, nops(R['variables']), R_Union);
    
    out := convertRankList(rankTable, R);
    
    return out;
    
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
          dim :: list(nonnegint),
          lrs_local :: TRDlrs,
          j :: nonnegint,
          i :: posint,
          lcs :: TRDlcs,
          lrs2 :: TRDlrs,
          rank :: table,
          tab :: table,
          perm :: list(posint);
          
    # Get the number of variables
    nVars := nops(R['variables']) - nParams;

    # Get the dimension of each regular system
    dim := map(getDimension, map(RC_CST:-RepresentingChain, lrs, R), nVars, R);

    # arrange in order of least dimension to greatest
    perm := sort(dim, 'output' = 'permutation');
    lrs_local := lrs[perm];

    # Convert each regular system to a constructible set
    lcs := map(RC:-TRDconstructible_set, map(`[]`, lrs_local),R);

    # Compute differences
    for i from 1 to nops(lcs) - 1 do
        # Compute difference with all constructible sets of higher dimension
        for j from i+1 to nops(lcs) do
            lcs[i] := RC:-TRDAlgebraicDifference(lcs[i], lcs[j], R);
        end do;
    end do;

    # Convert all constructible sets to regular systems
    lrs2 := map(x->op(RC_CST:-RepresentingRegularSystems(x, R)), lcs);

    # Construct the table with rank and regular systems (pre-projection)
    tab := constructTable(lrs2, nParams, R);

    # Initialize rank table
    rank := table([seq(i = {}, i = 0..nVars)]);

    # Project variables onto parameter space
    for i from 0 to nVars do
        if nops(tab[i]) > 0 then
            rank[i] := [linearProjection(convert(tab[i], list), nParams, R)];
        end if;
    end do;

    return rank;

end proc;


# ----------------------------------------------------------------------- #
# constructTable                                                          #
#                                                                         #
# Constructs a table such that table[r] corresponds to all regular        #
# systems of rank r from the input list.                                  #
#                                                                         #
# INPUT                                                                   #
#   lrs ....... a list of regular_systems                                 #
#   nParams ... number of parameters                                      #
#   R ......... polynomial ring                                           #
#                                                                         #
# OUTPUT                                                                  #
#   Returns a table where table[r] contains a set of regular_systems      #
#   corresponding to a rank of r.                                         #
# ----------------------------------------------------------------------- #
constructTable := proc(lrs::TRDlrs, nParams::posint, R::TRDring, $)

    local nVars :: posint,
          tab :: table,
          rs :: TRDrs,
          d :: nonnegint,
          i :: posint;

    # number of variables
    nVars := nops(R['variables']) - nParams;

    # initialize table
    tab := table();
    for i from 0 to nVars do
        tab[i] := {};
    end do;

    # Compute dimension and add to table
    for rs in lrs do
        d := getDimension(RC_CST:-RepresentingChain(rs, R), nVars, R);
        
        tab[nVars - d] := tab[nVars - d] union {rs};
    end do;

    for i from 1 to nVars do
        tab[i] := [op(tab[i])];
    end do;
    
    return tab;

end proc;


# ----------------------------------------------------------------------- #
# convertRankList                                                         #
#                                                                         #
# Converts a table where table[r] = CS where r is the computed rank and   #
# CS is a constructible set to a list with elements of the form [r, CS].  #
#                                                                         #
# INPUT                                                                   #
#    rankTable ... Table to convert                                       #
#    R ........... Polynomial ring                                        #
#                                                                         #
# OUTPUT                                                                  #
#   List with elements indicating the rank                                #
# ----------------------------------------------------------------------- #
convertRankList := proc(rankTable::table, R::TRDring, $)

    local idx :: list(nonnegint),
          i :: posint,
          outList :: list([nonnegint, TRDcs]);

    idx := sort([indices(rankTable,'nolist')]);

    outList := [];
    for i from 1 to nops(idx) do

        if nops(rankTable[idx[i]]) > 0 then
            outList := [op(outList), [idx[i], ListUnion(rankTable[idx[i]],R)]];
        end if;

    end do;

    return outList;

end proc;


# ----------------------------------------------------------------------- #
# linearProjection                                                        #
#                                                                         #
# Projects a list of regular systems onto the parameter space where all   #
# entries of lrs are linear in the unknown variables.                     #
#                                                                         #
# INPUT                                                                   #
#   lrs ....... a list of regular systems                                 #
#   nParams ... number of parameters                                      #
#   R ......... polynomial ring                                           #
#                                                                         #
# OUTPUT                                                                  #
#   A constructible set without the equations containing the unknowns.    #
# ----------------------------------------------------------------------- #
linearProjection := proc(lrs::TRDlrs, nParams::posint, R::TRDring, $)

    local nVars :: posint,
          R_params :: TRDring,
          ieqs,
          rc :: TRDrc,
          rs :: TRDrs,
          lrsUnder :: TRDlrs,
          rsUnder :: TRDrs,
          CS :: TRDcs;

    nVars := nops(R['variables']) - nParams;

    # new polynomial ring without variables
    R_params := RC:-PolynomialRing(R['variables'][nVars+1..-1]);

    # Under for a regular system
    lrsUnder := [];
    for rs in lrs do

        rc := RC_CST:-RepresentingChain(rs, R);

        ieqs := RC_CST:-RepresentingInequations(rs, R);

        rc := RC_CT:-Under(R['variables'][nVars], rc, R);

        rsUnder := RC_CST:-RegularSystem(rc, ieqs, R_params);

        lrsUnder := [op(lrsUnder), rsUnder];

    end do;
    
    lrsUnder := lrsUnder;

    CS := RC_CST:-ConstructibleSet(lrsUnder, R_params);

    return CS;

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