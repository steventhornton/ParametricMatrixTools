# ======================================================================= #
# ======================================================================= #
#                                                                         #
# comprehensive_rank.mpl                                                  #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Sept. 30/2017                                               #
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
        SuggestVariableOrderWithParameters,
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

    local n :: posint,
          m :: posint,
          y,
          lp,
          R_Union,
          CS_Int,
          lrs,
          rankTable,
          d,
          cs_singular,
          cs_nonsingular,
          out,
          i;

    # Get size of input matrix
    n, m := op(1,A);

    # Temporary variables for matrix
    y := seq('x'[i], i=1..m);
    
    # Get linear equations from matrix
    lp := convert(A.Vector([y]), list);

    # Join the parameters with the new variables.
    # Use suggest variable order for improved performance
    R_Union := RC:-PolynomialRing(SuggestVariableOrderWithParameters(lp, R['variables']));

    # Pre-condition for case where A is singular
    d := LA:-Determinant(A);
    
    # Parameter values where A is singular
    cs_singular := add_list_polys_to_cs([d], cs, R);
    
    # If the matrix is non-singular for all parameter values
    # return the input conditions with full rank
    if RC:-TRDAlgebraicIsEmpty(cs_singular, R) then
        return [[n, cs]];
    end if;
    
    # Parameter values where A is non-singular
    cs_nonsingular := RC:-TRDAlgebraicDifference(cs, cs_singular, R);
    
    # Add the equations from the matrix
    CS_Int := add_list_polys_to_cs(lp, cs_singular, R_Union);

    # Convert to list of regular systems
    lrs := RC:-TRDregular_systems(CS_Int, R_Union);

    # Call method to get rank
    rankTable := getRank(lrs, nops(R['variables']), R_Union);
    
    out := convertRankList(rankTable, R);
    
    if not RC:-TRDAlgebraicIsEmpty(cs_nonsingular, R) then
        out := [op(out), [n, cs_nonsingular]];
    end if;
    
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
          CS_tmp :: TRDcs,
          rank :: table,
          lrs2 :: TRDlrs,
          tab :: table;

    # Get the number of variables
    nVars := nops(R['variables']) - nParams;

    # Get the dimension of each regular system
    dim := map(getDimension, lrs, nVars, R);
    
    # arrange in order of least dimension to greatest
    lrs_local := [];
    for j from 0 to nVars do
        for i from 1 to nops(dim) do
            if evalb(dim[i] = j) then
                lrs_local := [op(lrs_local), lrs[i]];
            end if;
        end do;
    end do;

    lcs := map(RC:-TRDconstructible_set, map(`[]`, lrs_local),R);

    # Compute differences
    lrs2 := [];
    for i from 1 to nops(lrs_local) - 1 do

        # Create temporary constructible set
        CS_tmp := RC:-TRDconstructible_set([lrs_local[i]],R);

        # Compute difference with all regular systems later in the list
        for j from i+1 to nops(lrs_local) do
            CS_tmp := RC:-TRDAlgebraicDifference(CS_tmp, lcs[j], R);
        end do;

        # append to list of regular systems
         lrs2 := [op(lrs2), op(RC_CST:-RepresentingRegularSystems(CS_tmp, R))];

    end do;
    lrs2 := [op(lrs2), lrs_local[-1]];

    # Construct the table with rank and regular systems (pre-projection)
    tab := constructTable(lrs2, nParams, R);
    
    # initialize rank table
    rank := table();
    for i from 0 to nVars do
        rank[i] := {};
    end do;

    # project variables onto parameter space
    for i from 0 to nVars do
        if nops(tab[i]) > 0 then
            rank[i] := [linearProjection(convert(tab[i],list), nParams, R)];
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
        d := getDimension(rs, nVars, R);
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
# SuggestVariableOrderWithParameters                                      #
#                                                                         #
# Gives a variable ordering suggestion ensuring the parameters are less   #
# than any other variables.                                               #
#                                                                         #
# INPUT                                                                   #
#   lp ....... List of polynomials                                        #
#   params ... List of parameters                                         #
#                                                                         #
# OUTPUT                                                                  #
#   A list of variables.                                                  #
# ----------------------------------------------------------------------- #
SuggestVariableOrderWithParameters := proc(lp::list(polynom), params::list(name), $)
    
    local y::list(name), v::name;
    
    # Call standard suggest variable order
    y := RC:-SuggestVariableOrder(lp, params);
    
    # Make sure parameters are least variables
    for v in params do
        y := remove(x -> evalb(x=v), y);
    end do;
    return [op(LT:-Reverse(y)), op(params)];
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
# Computed the dimension (i.e. number of equations containing the         #
# unknowns) in a regular system.                                          #
#                                                                         #
# INPUT                                                                   #
#   rs ...... Regular system                                              #
#   nVars ... Number of variables                                         #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#   Returns a positive integer corresponding to the dimension of the      #
#   input rs.                                                             #
# ----------------------------------------------------------------------- #
getDimension := proc(rs::TRDrs, nVars::nonnegint, R::TRDring, $)

    local x :: set(name),
          rc :: TRDrc,
          eqns :: list(polynom),
          eqn :: polynom,
          i :: nonnegint;

    # Get the variables
    x := {op((R['variables'])[1..nVars])};

    # Get the equations of the regular system
    rc := RC_CST:-RepresentingChain(rs, R);
    eqns := RC:-Equations(rc, R);
    
    i := nVars;
    
    for eqn in eqns do
        if not evalb(nops(indets(eqn) intersect x) = 0) and not evalb(eqn = 0) then
            i := i - 1;
        end if;
    end do;
    
    return i;
    
end proc;

end module;