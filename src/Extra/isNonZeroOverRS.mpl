# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isNonZeroOverRS.mpl                                                     #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Determine if a polynomial vanishes nowhere in a regular system. That    #
# is, for a polynomial p, return true if                                  #
#       V(p) intersect {W(T)\V(h)}                                        #
# is the empty set and false otherwise. T is a regular chain and h is a   #
# polynomial representing the inequation constraints of the system where  #
# the input regular system represents all points in W(T)\V(h).            #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #                                                     #
#   rs ... Regular system                                                 #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   True if p vanishes nowhere in rs, false otherwise.                    #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([x, y]):                                        #
#   > rc := Chain([(x+y)^2], Empty(R), R):                                #
#   > rs := RegularSystem(rc, [y], R):                                    #
#   > p := y^2 + x + y:                                                   #
#   > isNonZeroOverRS(p, rs, R):                                          #
#         true                                                            #
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
isNonZeroOverRS := module()

    export ModuleApply;
    
    local
        init,
        implementation,
        isNonZeroOverRS_inequations,
        intersectRC,
        intersectCS;
    
    ModuleApply := proc(p::polynom, rs::TRDrs, R::TRDring, $) :: truefalse;
        return init(p, rs, R);
    end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ----------------------------------------------------------------------- #
# init                                                                    #
#                                                                         #
# Checks the types of the input and calls the implementation if all input #
# values pass checks.                                                     #
#                                                                         #
# INPUT                                                                   #
#   in_p ... Polynomial                                                   #
#   rs ..... Regular system                                               #
#   R ...... Polynomial ring                                              #
#                                                                         #
# OUTPUT                                                                  #
#   Same as isZeroOverRS                                                  #
# ----------------------------------------------------------------------- #
init := proc(in_p::polynom, rs::TRDrs, R::TRDring, $) :: truefalse;
    
    local p::polynom;
    
    # Ensure p is a polynomial in R
    if not RC:-TRDis_poly(in_p, R) then
        error "invalid polynomial: %1", in_p;
    else
        p := RC:-TRDmodularize_coefficients(in_p, R);
    end if;
    
    return implementation(p, rs, R);

end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Determine if a polynomial is non-zero everywhere over a regular system. #
#                                                                         #
# INPUT                                                                   #
#   in_p ... Polynomial                                                   #
#   rs ..... Regular system                                               #
#   R ...... Polynomial ring                                              #
#                                                                         #
# OUTPUT                                                                  #
#   - true if in_p is non-zero everywhere modulo rs                       #
#   - false otherwise                                                     #
# ----------------------------------------------------------------------- #
implementation := proc(in_p::polynom, rs::TRDrs, R::TRDring, $) :: truefalse;

    local p::polynom, 
          rc::TRDrc,
          h::list(polynom);
    
    # Reduce p modulo the saturated ideal of rc
    rc := RC:-TRDregular_chain_rs(rs, R);
    p := RC:-SparsePseudoRemainder(in_p, rc, R);
    
    # Case where p is constant and non-zero
    if RC:-TRDis_constant(p, R) then
        return evalb(p <> 0);
    end if;
    
    # Check if the input polynomial divides the inequations of rs
    if isNonZeroOverRS_inequations(p, rs, R) then
        return true;
    end if;
    
    h := RC_CST:-RepresentingInequations(rs, R);
    
    if nops(h) = 0 then
        return intersectRC(p, rs, R);
    else
        return intersectCS(p, rs, R);
    end if;

end proc;


# ----------------------------------------------------------------------- #
# isNonZeroOverRS_inequations                                             #
#                                                                         #
# Check if a polynomial is non-zero over the inequations of a regular     #
# system.                                                                 #
#                                                                         #
# INPUT                                                                   #
#   in_p ... Polynomial                                                   #
#   rs ..... Regular system                                               #
#   R ...... Polynomial ring                                              #
#                                                                         #
# OUTPUT                                                                  #
#   - True if in_p and the inequation are not co-prime                    #
#   - False otherwise                                                     #
#                                                                         #
# ASSUMPTIONS                                                             #
#   p is assumed to be reduced modulo rc and expanded                     #
# ----------------------------------------------------------------------- #
isNonZeroOverRS_inequations := proc(in_p::polynom, rs::TRDrs, R::TRDring, $) :: truefalse;
    
    local h1::{list(polynom), set(polynom)},
          h2::{list(polynom), set(polynom)},
          h::polynom,
          rc::TRDrc,
          p::{polynom, list(polynom)},
          g::polynom;
    
    # Get the regular chain from rs
    rc := RC:-TRDregular_chain_rs(rs, R);
    
    # Get the inequations from rs
    h1 := RC:-TRDinequations_rs(rs, R);
    h2 := RC:-Inequations(rc, R); 
    
    # Multiply all inequations together
    h := RC:-TRDlist_mul_polys([1, op(h1), op(h2)], R);
    
    # Make p square-free
    p := RC:-TRDGcdFreeFactorization(in_p, R);
    p := RC:-TRDlist_mul_polys(p, R);
    
    # Compute the gcd of the inequation and the polynomial
    g := gcd(h, p);
    
    return evalb(expand(g-p) = 0)

end proc;


# ----------------------------------------------------------------------- #
# intersectRC                                                             #
#                                                                         #
# [Description]                                                           #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   rs ... Regular system                                                 #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#                                                                         #
# ASSUMPTIONS                                                             #
#   p is assumed to be reduced modulo rc and expanded                     #
# ----------------------------------------------------------------------- #
intersectRC := proc(p::polynom, rs::TRDrs, R::TRDring, $) :: truefalse;
    
    local rc::TRDrc,
          lrc::list(TRDrc);
    
    # Get the regular chain from rs
    rc := RC:-TRDregular_chain_rs(rs, R);
    
    # Need a faster way of doing this...
    # Remove extension phase of intersect?
    lrc := RC:-Intersect(p, rc, R);
    
    return evalb(nops(lrc) = 0);

end proc;


# ----------------------------------------------------------------------- #
# intersectCS                                                             #
#                                                                         #
# [Description]                                                           #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   rs ... Regular system                                                 #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#                                                                         #
# ASSUMPTIONS                                                             #
#   p is assumed to be reduced modulo rc and expanded                     #
# ----------------------------------------------------------------------- #
intersectCS := proc(p::polynom, rs::TRDrs, R::TRDring, $) :: truefalse;
    
    local cs1::TRDcs,
          cs2::TRDcs,
          csI::TRDcs;
    
    # Convert rs into a constructible set
    cs1 := RC_CST:-ConstructibleSet([rs], R);
    
    # Convert p = 0 into a constructible set
    cs2 := RC_CST:-GeneralConstruct([p], [], R);
    
    # Compute intersection
    csI := RC_CST:-Intersection(cs1, cs2, R);
    
    # True if intersection is empty, occurs when p <> 0 in
    # rs, we found a contradiction in our assumption in cs2 that p = 0
    return RC_CST:-IsEmpty(csI, R);

end proc;

end module;