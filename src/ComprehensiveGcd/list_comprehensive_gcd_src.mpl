# ======================================================================= #
# ======================================================================= #
#                                                                         #
# list_comprehensive_gcd_src.mpl                                          #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 24/2017                                                #
#                                                                         #
# Computes the gcd of a list of parametric univariate polynomials in the  #
# sense of Lazard. Computation is done over a constructible set.          #
#                                                                         #
# INPUT                                                                   #
#   lp ... List of polynomials                                            #
#   v .... Variable                                                       #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A sequence gcdList, cs_zero where gcdList is a list of pairs of the   #
#   form:                                                                 #
#       [g_i, cs_i]                                                       #
#   where g_i is the gcd of the polynomials in lp for all values in the   #
#   zero set of cs_i. cs_zero is the constructible set where both more    #
#   than one polynomial in lp vanishes for all values in its zero set.    #
#   The set {cs_1, cs_2, ..., cs_zero} forms a partition of the input     #
#   constructible set.                                                    #
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
list_comprehensive_gcd_src := module()

    export ModuleApply;
    
    local
        implementation,
        sortListPureLex,
        listGcd;
    
    ModuleApply := proc(lp::depends(list(polyInRing(R))), v::name, cs::TRDcs, R::TRDring, $) :: list([polynom, TRDcs]), TRDcs;
        return implementation(lp, v, cs, R);
    end proc:

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Compute the gcd of a list of polynomials over a constructible set.      #
#                                                                         #
# INPUT/OUTPUT                                                            #
#    Same as list_comprehensive_gcd_src                                   #
# ----------------------------------------------------------------------- #
implementation := proc(lp_in::depends(list(polyInRing(R))), v::name, cs::TRDcs, R::TRDring, $)
    
    local lp :: list(polynom),
          gList :: list([polynom, TRDcs]),
          cs_zero :: TRDcs,
          gStack :: Stack,
          pStack :: Stack,
          nextStack :: Stack,
          i :: posint,
          p :: polynom,
          g :: polynom,
          es :: TRDcs,
          gNew :: list([polynom, TRDcs]),
          cs_zero_new :: TRDcs;
    
    # If all polynomials in lp only contain v as an indeterminant, call 
    # non-parametric method
    if nops(indets(lp_in)) = 0 then
        return [listGcd(lp_in), cs];
    end if;
    if nops(indets(lp_in)) = 1 and v in indets(lp_in) then
        return [listGcd(lp_in), cs];
    end if;

    # Order the elements of lp in decreasing order
    lp := sortListPureLex(lp_in, R);

    # Compute gcd of first two elements of lp
    gList, cs_zero := ComprehensiveGcd(lp[1], lp[2], v, cs, R, 'cofactors'=false, 'outputType'='CS');
    
    # If there is only 2 elements, return
    if nops(lp) = 2 then
        return gList, cs_zero;
    end if;
    
    # --------------------------- #
    # Compute the rest of the gcd #
    # --------------------------- #
    
    gStack := SimpleStack();        # [polynom, TRDcs]
    pStack := SimpleStack();        # polynom
    nextStack := SimpleStack();     # [polynom, TRDcs]
    
    # Fill the gcd stack
    for i to nops(gList) do
        gStack:-push(gList[i]);
    end do;
    
    # Add polynomials that have not been considered in the gcd yet to the stack
    for i from 3 to nops(lp) do
        pStack:-push(lp[i]);
    end do;
    
    # Add each polynomial to gcd computation one at a time
    while not pStack:-empty() do
        
        # Get the next polynomial
        p := pStack:-pop();
        
        # Compute gcd for each branch of computation
        while not gStack:-empty() do
        
            g, es := op(gStack:-pop());
            
            gNew, cs_zero_new := ComprehensiveGcd(p, g, v, es, R, 'cofactors'=false, 'outputType'='CS');
            
            # Update cs_zero
            cs_zero := RC_CST:-Union(cs_zero, cs_zero_new, R);
            
            # Add each element of gNew to nextStack
            for i to nops(gNew) do
                nextStack:-push(gNew[i]);
            end do;
        
        end do;
        
        # gStack is now empty.
        # move everything from nextStack to gStack
        while not nextStack:-empty() do
            gStack:-push(nextStack:-pop());
        end do;
        
    end do;
    
    # Convert gStack to list and return
    return [seq(gStack:-pop(), i=1..gStack:-depth())], cs_zero;

end proc;


# ----------------------------------------------------------------------- #
# sortListPureLex                                                         #
#                                                                         #
# Method to sort (decreasing order) the elements of a list by their       #
# initials using pure lexicographical ordering.                           #
#                                                                         #
# INPUT                                                                   #
#    lp_in ... A list of polynomials                                      #
#    R ....... Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#    A list with the elements sorted                                      #
# ----------------------------------------------------------------------- #
sortListPureLex := proc(lp_in::depends(list(polyInRing(R))), R::TRDring, $) :: list(polynom);
    return LT:-Reverse(sort(lp_in, proc(x, y) options operator, arrow; RC:-TRDstrictly_less_ritt(x, y, R) end proc, 'output'='sorted'));
end proc;


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# External Files
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveGcd/listGcd.mpl>

end module;