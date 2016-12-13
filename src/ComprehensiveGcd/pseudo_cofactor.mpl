# ======================================================================= #
# ======================================================================= #
#                                                                         #
# pseudo_cofactor.mpl                                                     #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 12/2016                                                #
#                                                                         #
# Computes the polynomial p/g where g is known to be a factor of p over   #
# a regular chain.                                                        #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   g .... A factor (gcd) of p                                            #
#   v .... Variable                                                       #
#   rc ... Regular chain                                                  #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   The polynomial p/g mod rc.                                            #
#                                                                         #
# ASSUMPTIONS                                                             #
#   g is never 0                                                          #
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
pseudo_cofactor := proc(p_in::depends(polyInRing(R)), g_in::depends(polyInRing(R)), v::name, rc::TRDrc, R::TRDring, $)
    
    local r :: polynom,
          m :: polynom,
          q :: ratpoly,
          p :: polynom,
          g :: polynom;
    
    p := RC:-SparsePseudoRemainder(p_in, rc, R);
    g := RC:-SparsePseudoRemainder(g_in, rc, R);
    
    ASSERT(not isZeroOverRS(g, RC_CST:-RegularSystem(rc, [], R), R), "g must not be zero");
    
    r := sprem(p, g, v, 'm','q');
    
    ASSERT(isZeroOverRS(r, RC_CST:-RegularSystem(rc, R), R), "Remainder must be zero");
    
    q := RC:-SparsePseudoRemainder(q, rc, R);
    q := normal(q/m);
    
    return q;

end proc: