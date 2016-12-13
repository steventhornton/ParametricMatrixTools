# ======================================================================= #
# ======================================================================= #
#                                                                         #
# compute_cofactors.mpl                                                   #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 13/2016                                                #
#                                                                         #
# Computes the co-factors of two polynomials. That is, the rational       #
# polynomials cof_p1 = p1/g and cof_p2 = p2/g where g = gcd(p1, p2).      #
#                                                                         #
# INPUT                                                                   #
#   p1 ... Polynomial                                                     #
#   p2 ... Polynomial                                                     #
#   g .... The gcd of p1 and p2 such that it is the gcd of p1 and p2      #
#          for all values in the zero set of rs.                          #
#   v .... Variable                                                       #
#   rs ... Regular system                                                 #
#   R .... Polynomial Ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list of the form                                                    #
#       [g, cof_a, cof_b, rs]                                             #
#   g     : Polynomial                                                    #
#   cof_a : Cofactor of p1                                                #
#   cof_b : Cofactor of p2                                                #
#   rs    : Regular system                                                #
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
compute_cofactors_rs := proc(p1::depends(polyInRing(R)), p2::depends(polyInRing(R)), g::depends(polyInRing(R)), v::name, rs::TRDrs, R::TRDring, $)
    
    local rc :: TRDrc,
          q1 :: ratpoly,
          q2 :: ratpoly;
    
    ASSERT(not isZeroOverRS(g, rs, R), "g must not be zero");
    
    rc := RC_CST:-RepresentingChain(rs, R);
    
    q1 := pseudo_cofactor(p1, g, v, rc, R);
    q2 := pseudo_cofactor(p2, g, v, rc, R);
    
    ASSERT(isZeroOverRS(numer(g*q1 - p1), rs, R), "Cofactor 1 not computed correctly");
    ASSERT(isZeroOverRS(numer(g*q2 - p2), rs, R), "Cofactor 2 not computed correctly");
    
    return [q1, q2];
    
end proc;