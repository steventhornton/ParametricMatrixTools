# ======================================================================= #
# ======================================================================= #
#                                                                         #
# comprehensive_gcd_src.mpl                                               #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 9/2016                                                 #
#                                                                         #
# Computes the gcd of two parametric univariate polynomials in the sense  #
# of Lazard. Computation is done over a constructible set.                #
#                                                                         #
# INPUT                                                                   #
#   p1 ... Polynomial                                                     #
#   p2 ... Polynomial                                                     #
#   v .... Variable                                                       #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A sequence gcdList, cs_zero where gcd_list is a list of pairs of the  #
#   form:                                                                 #
#       [g_i, cs_i]                                                       #
#   where g_i is the gcd of p1 and p2 for all values in the zero set of   #
#   cs_i. cs_zero is the constructible set where both p1 and p2 vanish    #
#   for all values in its zero set. The set {cs_1, cs_2, ..., cs_zero}    #
#   forms a partition of the input constructible set.                     #
#                                                                         #
# REFERENCES                                                              #
#                                                                         #
# TO DO                                                                   #
#   - Move cs_zero to one call                                            #
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
comprehensive_gcd_src := module()

    export ModuleApply;

    local implementation,
          gcd_for_zero_resultant_cs,
          gcd_for_zero_resultant_rs,
          hasZeroPoly,
          gcd_for_constants,
          missingMainVar,
          get_gcd_no_subresultant,
          gcd_by_subresultant_non_vanishing_initals;

    ModuleApply := proc(p1::depends(polyInRing(R)), 
                        p2::depends(polyInRing(R)), 
                        v::name, 
                        cs::TRDcs,
                        R::TRDring, $)::list([polynom, TRDcs]), TRDcs;
        
        local result, cs_zero;
        
        # Compute the gcd
        result, cs_zero := implementation(p1, p2, v, cs, R);
        
        # Remove empty elements and entries with empty constructible sets
        result := remove(x -> nops(x) = 0, result);
        result := remove(x -> RC:-TRDis_empty_constructible_set(x[2],R), result);
        
        return result, cs_zero;
        
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Compute the gcd of two polynomials over a constructible set.            #
#                                                                         #
# INPUT/OUTPUT                                                            #
#    Same as comprehensive_gcd_src                                        #
# ----------------------------------------------------------------------- #
implementation := proc(p1, p2, v, cs, R)
    
    local cs_nz :: TRDcs,
          cs_z :: TRDcs,
          cs_zero :: TRDcs,
          result :: list({[polynom, TRDcs], []}),
          result_tmp :: list({[polynom, TRDcs], []}),
          r :: polynom,
          src :: TRDsrc; 
    
    ASSERT(RC:-TRDis_expanded(p1));
    ASSERT(RC:-TRDis_expanded(p2));
    
    # Check if cs is an empty constructible set
    if RC:-TRDis_empty_constructible_set(cs, R) then
        return [[]], RC:-TRDempty_constructible_set();
    end if;
    
    # Check if one of p1 or p2 contain no variables
    if RC:-TRDis_constant(p1, R) or RC:-TRDis_constant(p2, R) then
        return gcd_for_constants(p1, p2, v, cs, R);
    end if;
    
    # Check if one of p1 or p2 vanish at all points in the zero set of cs
    if isZeroOverCS(p1, cs, R) or isZeroOverCS(p2, cs, R) then
        return hasZeroPoly(p1, p2, cs, R);
    end if;
    
    # Check if either p1 or p2 is constant w.r.t. v
    if RC:-TRDuniv_degree(p1, v) = 0 or RC:-TRDuniv_degree(p2, v) = 0  then
        return missingMainVar(p1, p2, v, cs, R);
    end if;
    
    # Compute the subresultant chain and get the resultant
    src := RC_CT:-SubresultantChain(p1, p2, v, R);
    r := RC_CT:-SubresultantOfIndex(0, src, R);
    
    # Split into cases 2 cases cs_nz and cs_z that form a partition of cs:
    #   cs_nz: For all points in the zero set of cs, r is non-zero
    #    cs_z: r vanishes at all points in the zero set of cs
    cs_nz, cs_z := TRDdifference_intersect_cs_p(cs, r, R);
    
    # res(p1, p2) <> 0 <=> gcd(p1, p2) = 1
    result := [[1, cs_nz]];
    
    # Call method for the case where res(p1, p2) vanishes at all points in 
    # the zero set of cs.
    result_tmp, cs_zero := gcd_for_zero_resultant_cs(p1, p2, src, v, cs_z, R);
    result := [op(result), op(result_tmp)];
    
    return result, cs_zero;
    
end proc;


# ----------------------------------------------------------------------- #
# gcd_for_zero_resultant_cs                                               #
#                                                                         #
# Compute the gcd of two polynomials given a constructible set such that  #
# res(p1, p2) = 0 at all points in the zero set of the input              #
# constructible set.                                                      #
#                                                                         #
# INPUT/OUTPUT                                                            #
#    Same as comprehensive_gcd_src                                        #
# ----------------------------------------------------------------------- #
gcd_for_zero_resultant_cs := proc(p1, p2, src, v, cs, R, $)
    
    local lrs::TRDlrs,
          rs::TRDrs, 
          result :: list({[polynom, TRDcs], []}),
          cs_zero::TRDcs, 
          result_tmp :: list({[polynom, TRDcs], []}),
          cs_tmp::TRDcs, 
          cs_disjoint::TRDcs;
    
    ASSERT(RC:-TRDis_expanded(p1));
    ASSERT(RC:-TRDis_expanded(p2));
    
    # Initialize the output
    result := [];
    cs_zero := RC:-TRDempty_constructible_set();
    
    # Make cs disjoint (may not be neccessary)
    cs_disjoint := RC:-TRDconstructible_set_make_pairwise_disjoint(cs, R);
    
    # Call gcd_for_zero_resultant_rs on each regular system in cs
    lrs := RC:-TRDregular_systems(cs_disjoint, R);
    for rs in lrs do
        result_tmp, cs_tmp := gcd_for_zero_resultant_rs(p1, p2, src, v, rs, R);
        result := [op(result), op(result_tmp)];
        cs_zero := RC_CST:-Union(cs_zero, cs_tmp, R);
    end do;
    
    return result, cs_zero;
    
end proc;


# ----------------------------------------------------------------------- #
# gcd_for_zero_resultant_rs                                               #
#                                                                         #
# Compute the gcd of two polynomials given a regular system such that     #
# res(p1, p2) = 0 at all points in the zero set of the input regular      #
# system.                                                                 #
#                                                                         #
# INPUT                                                                   #
#   p1 .... Polynomial                                                    #
#   p2 .... Polynomial                                                    #
#   src ... Subresultant chain of p1 and p2                               #
#   v ..... Variable                                                      #
#   rs .... Regular system                                                #
#   R ..... Polynomial ring                                               #
#                                                                         #
# OUTPUT                                                                  #
#    Same as comprehensive_gcd_src                                        #
#                                                                         #
# ASSUMPTIONS                                                             #
#   deg(p1, v) > 0                                                        #
#   deg(p2, v) > 0                                                        #
#   res(p1, p2, v) = 0 in the zero set of rs                              #
# ----------------------------------------------------------------------- #
gcd_for_zero_resultant_rs := proc(p1::depends(polyInRing(R)), p2::depends(polyInRing(R)), src::TRDsrc, v::name, rs::TRDrs, R::TRDring, $)
    
    local p1_init::polynom,
          p2_init::polynom,
          cs_init::TRDcs,
          cs::TRDcs,
          cs_nz::TRDcs,
          cs_z::TRDcs,
          p1_tail::polynom,
          p2_tail::polynom,
          result::list({[polynom, TRDcs], []}),
          result_tmp::list({[polynom, TRDcs], []}),
          cs_zero::TRDcs,
          cs_zero_tmp::TRDcs;
    
    ASSERT(RC:-TRDis_expanded(p1));
    ASSERT(RC:-TRDis_expanded(p2));
    
    cs := RC:-TRDconstructible_set([rs], R);
    
    # Split cs into two constructible sets:
    #   cs_nz: For all points in the zero set of cs, the initials of p1 and 
    #          p2 don't simultaneously vanish
    #    cs_z: init(p1) and init(p2) vanish at all points in the zero set 
    #          of cs_z
    p1_init := RC:-TRDuniv_lcoeff(p1, v);
    p2_init := RC:-TRDuniv_lcoeff(p2, v);
    cs_init := RC_CST:-GeneralConstruct([p1_init, p2_init], [], R);
    cs_nz, cs_z:= RC:-TRDdifference_intersect_cs_cs(cs, cs_init, R);
    
    # Recursive call on Tail(p1), Tail(p2)
    p1_tail := RC:-TRDuniv_tail(p1, v);
    p2_tail := RC:-TRDuniv_tail(p2, v);
    result, cs_zero := comprehensive_gcd_src(p1_tail, p2_tail, v, cs_z, R);
    
    # Get the resultant in the case where the initials of p1 and p2 don't 
    # simultaneously vanish
    result_tmp, cs_zero_tmp := gcd_by_subresultant_non_vanishing_initals(p1, p2, src, v, cs_nz, R);
    result := [op(result), op(result_tmp)];
    
    cs_zero := RC_CST:-Union(cs_zero, cs_zero_tmp, R);
    
    return result, cs_zero;
    
end proc;


# ----------------------------------------------------------------------- #
# gcd_by_subresultant_non_vanishing_initals                               #
#                                                                         #
# Compute the gcd by finding the first subresultant whoes initial does    #
# not vanish on cs given that the resultant of p1 and p2 vanishes         #
# everywhere in the zero set of cs_in, and the initials of p1 and p2      #
# don't simultaneously vanish.                                            #
#                                                                         #
# INPUT/OUTPUT                                                            #
#    Same as comprehensive_gcd_src                                        #
# ----------------------------------------------------------------------- #
gcd_by_subresultant_non_vanishing_initals := proc(p1::depends(polyInRing(R)), p2::depends(polyInRing(R)), src::TRDsrc, v::name, cs_in::TRDcs, R::TRDring, $)
    
    local cs::TRDcs,
          g::polynom,
          i::posint,
          init_g::polynom,
          cs_nz::TRDcs,
          result::list({[polynom, TRDcs], []}),
          cs_zero::TRDcs,
          result_tmp::list({[polynom, TRDcs], []});
    
    cs := RC:-TRDrename_constructible_set(copy(cs_in));
    
    result := [];
    
    # Proceeding bottom up with src, get the first subresultant whos 
    # initial does not vanish in the zero set of cs; it is the gcd. If no 
    # such subresultant exists, either p1 or p2 is the gcd.
    for i to nops(src['subresultant_chain_vector']) -2 - 1 while not RC:-TRDis_empty_constructible_set(cs, R) do
        
        g := RC_CT:-SubresultantOfIndex(i, src, R);
        
        # The computation is split into two cases: 
        #   1. The lcoeff(g, v) does not vanish anywhere in the zero set of
        #      cs_nz, and thus g is the gcd of p1 and p2.
        #   2. The lcoeff(g, v) vanishes everywhere in the zero set of cs, 
        #      computation is continued on the next subresultant.
        init_g := RC:-TRDuniv_lcoeff(g, v);
        cs_nz, cs := TRDdifference_intersect_cs_p(cs, init_g, R);
        
        result := [op(result), [g, cs_nz]];
        
    end do;
    
    # If cs is non-empty, either p1 or p2 is the gcd
    result_tmp, cs_zero := get_gcd_no_subresultant(p1, p2, v, cs, R);
    result := [op(result), op(result_tmp)];
    return result, cs_zero;

end proc:


# ----------------------------------------------------------------------- #
# get_gcd_no_subresultant                                                 #
#                                                                         #
# Get the gcd of p1 and p2 when none of the subresultants of p1 and p2    #
# have an initial that does not vanish in the zero set of a constructible #
# set. In this case, either p1 or p2 is the gcd.                          #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as comprehensive_gcd_src                                         #
# ----------------------------------------------------------------------- #
get_gcd_no_subresultant := proc(p1, p2, v, cs, R)
    
    local g::polynom,
          cs_p1 :: TRDcs,
          cs_p2 :: TRDcs,
          cs_p1_nz :: TRDcs,
          cs_p1_z :: TRDcs,
          cs_p2_nz :: TRDcs,
          cs_z_z :: TRDcs,
          cs_z_nz :: TRDcs,
          cs_nz_z :: TRDcs,
          cs_nz_nz :: TRDcs;
    
    if RC:-TRDis_empty_constructible_set(cs, R) then
        return [[]], RC:-TRDempty_constructible_set();
    end if;
    
    # Need to ensure neither p1 nor p2 are zero.
    if degree(p1, v) < degree(p2, v) then
        g := p1;
    else
        g := p2;
    end if;
    
    # Get cases where p1 and p2 are both zero
    cs_p1 := RC_CST:-GeneralConstruct([coeffs(p1, v)], [], R);
    cs_p2 := RC_CST:-GeneralConstruct([coeffs(p2, v)], [], R);
    
    cs_p1_nz, cs_p1_z := RC:-TRDdifference_intersect_cs_cs(cs, cs_p1, R);
    cs_p2_nz := RC_CST:-Difference(cs, cs_p2, R);
    
    cs_z_z, cs_z_nz := RC:-TRDdifference_intersect_cs_cs(cs_p1_z, cs_p2_nz, R);
    cs_nz_z, cs_nz_nz := RC:-TRDdifference_intersect_cs_cs(cs_p1_nz, cs_p2_nz, R);
    
    return [[g, cs_nz_nz], [p1, cs_nz_z], [p2, cs_z_nz]], cs_z_z;
    
end proc;


# ----------------------------------------------------------------------- #
# hasZeroPoly                                                             #
#                                                                         #
# Compute the gcd of p1 and p2 given that one of either p1 or p2 vanishes #
# everywhere in the zero set of a constructible set.                      #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as comprehensive_gcd_src                                         #
# ----------------------------------------------------------------------- #
hasZeroPoly := proc(p1::depends(polyInRing(R)), p2::depends(polyInRing(R)), cs::TRDcs, R::TRDring, $)
    
    local p1_z::truefalse, p2_z::truefalse;
    
    ASSERT(RC:-TRDis_expanded(p1));
    ASSERT(RC:-TRDis_expanded(p2));
    ASSERT(evalb(not RC:-TRDis_empty_constructible_set(cs, R)), "cs must not be empty.");
    
    p1_z := isZeroOverCS(p1, cs, R);
    p2_z := isZeroOverCS(p2, cs, R);
    
    if p1_z and p2_z then
        return [[]], cs;
    end if;
    
    if p1_z then
        return [[p2, cs]], RC:-TRDempty_constructible_set();
    else
        return [[p1, cs]], RC:-TRDempty_constructible_set();
    end if;
    
end proc;


# ----------------------------------------------------------------------- #
# gcd_for_constants                                                       #
#                                                                         #
# Compute the gcd of p1 and p2 when one of either p1 or p2 contains       #
# no indeterminants.                                                      #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as comprehensive_gcd_src                                         #
# ----------------------------------------------------------------------- #
gcd_for_constants := proc(p1::depends(polyInRing(R)), p2::depends(polyInRing(R)), v::name, cs::TRDcs, R::TRDring, $)
    
    local cs_p1 :: TRDcs,
          cs_nz :: TRDcs,
          cs_z :: TRDcs;
    
    ASSERT(RC:-TRDis_expanded(p1));
    ASSERT(RC:-TRDis_expanded(p2));
    ASSERT(evalb(not RC:-TRDis_empty_constructible_set(cs, R)), "cs must not be empty.");
    
    # Ensure p2 contains no variables
    if not RC:-TRDis_constant(p2, R) then
        #A recursive call is made to ensure p2 is constant, p1 may or may 
        # not be constant.
        return gcd_for_constants(p2, p1, v, cs, R);
    end if;
    
    # Case where both p1 and p2 are constant
    if RC:-TRDis_constant(p1, R) then
        if p1 = 0 and p2 = 0 then
            return [[]], cs;
        else
            return [[gcd(p1, p2), cs]], RC:-TRDempty_constructible_set();
        end if;
    end if;
    
    # Case where p2 = 0
    if p2 = 0 then
        # Case where both p1 and p2 are 0
        if isZeroOverCS(p1, cs, R) then
            return [[]], cs;
        end if;
        
        # Split for cases where p1=0 and p1<>0
        cs_p1 := RC_CST:-GeneralConstruct([coeffs(expand(p1), v)], [], R);
        cs_nz, cs_z := RC:-TRDdifference_intersect_cs_cs(cs, cs_p1, R);
        
        return [[p1, cs_nz]], cs_z;
        
    end if;
    
    # If p2 <> 0, call gcd
    return [[gcd(p1, p2), cs]], RC:-TRDempty_constructible_set();
    
end proc;


# ----------------------------------------------------------------------- #
# missingMainVar                                                          #
#                                                                         #
# Compute the gcd of p1 and p2 when one of either p1 or p2 does not       #
# contain the variable the gcd is being computed w.r.t.                   #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as comprehensive_gcd_src                                         #
# ----------------------------------------------------------------------- #
missingMainVar := proc(p1::depends(polyInRing(R)), p2::depends(polyInRing(R)), v::name, cs::TRDcs, R::TRDring, $)

    local cs_p1_coeff::TRDcs,
          cs_p2::TRDcs,
          cs1_z::TRDcs,
          cs1_nz::TRDcs,
          cs2_nz::TRDcs,
          cs_z_z::TRDcs,
          cs_z_nz::TRDcs,
          cs_nz_z::TRDcs,
          cs_nz_nz::TRDcs,
          cs_nz::TRDcs, 
          cs_z::TRDcs, 
          cs_p12::TRDcs;
    
    ASSERT(RC:-TRDis_expanded(p1));
    ASSERT(RC:-TRDis_expanded(p2));
    ASSERT(evalb(not RC:-TRDis_empty_constructible_set(cs, R)), "cs must not be empty.");
    
    # Ensure v >= mvar(p1) >= mvar(p2) and v > mvar(p2).
    if RC:-MainVariable(p2, R) = v then
        return missingMainVar(p2, p1, v, cs, R);
    end if;
    
    # If neither p1 nor p2 are polynomials in v, return splitting case
    if RC:-MainVariable(p1, R) <> v then
        cs_p12 := RC_CST:-GeneralConstruct([p1, p2], [], R);
        cs_nz, cs_z := RC:-TRDdifference_intersect_cs_cs(cs, cs_p12, R);
        return [[gcd(p1, p2), cs_nz]], cs_z;
    end if;
    
    # If mvar(p1) = v, split into cases where p1 and p2 <> 0
    cs_p1_coeff := RC_CST:-GeneralConstruct([coeffs(p1, v)], [], R);
    cs_p2 := RC_CST:-GeneralConstruct([p2], [], R);
    
    cs1_nz, cs1_z := RC:-TRDdifference_intersect_cs_cs(cs, cs_p1_coeff, R);
    
    cs2_nz := RC_CST:-Difference(cs, cs_p2, R);
    
    cs_z_z, cs_z_nz := RC:-TRDdifference_intersect_cs_cs(cs1_z, cs2_nz, R);
    cs_nz_z, cs_nz_nz := RC:-TRDdifference_intersect_cs_cs(cs1_nz, cs2_nz, R);
    
    return [[p2, cs_z_nz], [p1, cs_nz_z], [1, cs_nz_nz]], cs_z_z;
    
end proc;

end module;