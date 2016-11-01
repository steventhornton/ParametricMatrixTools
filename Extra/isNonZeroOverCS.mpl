# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isNonZeroOverCS.mpl                                                     #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Determine if a polynomial vanishes nowhere in a constructible set. That #
# is, for a polynomial p, return true if for each regular system [T, h]   #
# in a constructible set                                                  #
#       V(p) intersect {W(T)\V(h)}                                        #
# is the empty set. T is a regular chain and h is a polynomial            #
# representing the inequation constraints of the system where the input   #
# regular system represents all points in W(T)\V(h).                      #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #                                                     #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   True if p vanishes nowhere in cs, false otherwise.                    #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([x, y]):                                        #
#   > cs := GeneralConstruct([x+1], [], R):                               #
#   > p := x + y + 1:                                                     #
#   > isNonZeroOverCS(p, cs, R);                                          #
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
isNonZeroOverCS := module()

    export ModuleApply;

    local
        init,
        implementation;

    ModuleApply := proc(p::polynom, cs::TRDcs, R::TRDring, $) :: truefalse;
        init(p, cs, R);
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
#   cs ..... Constructible set                                            #
#   R ...... Polynomial ring                                              #
#                                                                         #
# OUTPUT                                                                  #
#   Same as isNonZeroOverCS                                               #
# ----------------------------------------------------------------------- #
init := proc(in_p::polynom, cs::TRDcs, R::TRDring, $) :: truefalse;
    
    local p::polynom;
    
    # Ensure p is a polynomial in R
    if not RC:-TRDis_poly(in_p, R) then
        error "invalid polynomial: %1", in_p;
    else
        p := RC:-TRDmodularize_coefficients(in_p, R);
    end if;
    
    return implementation(p, cs, R);
    
end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Determine if a polynomial is nonzero everywhere modulo the equations    #
# and inequations in a constructible set.                                 #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   Same as isNonZeroOverCS                                               #
# ----------------------------------------------------------------------- #
implementation := proc(p::polynom, cs::TRDcs, R::TRDring, $) :: truefalse;
    
    local rs::TRDrs;
    
    # Case where p is constant and non-zero
    if RC:-TRDis_constant(p, R) then
        return evalb(p <> 0);
    end if;
    
    # p must be non-zero zero over all representing regular systems in cs
    for rs in RC_CST:-RepresentingRegularSystems(cs, R) do
        if not isNonZeroOverRS(p, rs, R) then
            return false;
        end if;
    end do;
    
    return true;
    
end proc;

end module;