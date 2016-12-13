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
# Computes the co-factors of two polynomials. That is, the polynomials    #
# cof_p1 = p1/g and cof_p2 = p2/g where g = gcd(p1, p2).                  #
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
compute_cofactors_rs := proc(p1_in::depends(polyInRing(R)), p2_in::depends(polyInRing(R)), g::depends(polyInRing(R)), v::name, rs::TRDrs, R::TRDring, $)
    
    local rc :: TRDrc,
          q1 :: ratpoly,
          q2 :: ratpoly, qMonic, p1, p2, p1_lc, p2_lc, h1, h2;
    
    ASSERT(not isZeroOverRS(g, rs, R), "g must not be zero");
    
    rc := RC_CST:-RepresentingChain(rs, R);
    
    # p1 := RC:-SparsePseudoRemainder(p1_in, rc, R, 'h1');
    # p2 := RC:-SparsePseudoRemainder(p2_in, rc, R, 'h2');
    
    # p1_lc := lcoeff(RC:-SparsePseudoRemainder(p1, rc, R), v);
    # if isNonZeroOverRS(p1_lc, rs, R) then
    #     qMonic := normal(p1/p1_lc);
    #     if denom(qMonic) = 1 then
    #         p1 := qMonic;
    #     end if;
    # end if;
    # p2_lc := lcoeff(RC:-SparsePseudoRemainder(p2, rc, R), v);
    # print(p2_lc);
    # if isNonZeroOverRS(p2_lc, rs, R) then
    #     qMonic := normal(p2/p2_lc);
    #     if denom(qMonic) = 1 then
    #         p2 := qMonic;
    #     end if;
    # end if;
    # print(p2);
    
    q1 := pseudo_cofactor(p1_in, g, v, rc, R);
    q2 := pseudo_cofactor(p2_in, g, v, rc, R);
    
    # q1 := normal(q1/h1);
    # q2 := normal(q2/h2);
    
    # if isNonZeroOverRS(lcoeff(q1, v), rs, R) then
    #     qMonic := normal(q1/lcoeff(q1, v));
    #     if denom(qMonic) = 1 then
    #         q1 := qMonic;
    #     end if;
    # end if;
    # if isNonZeroOverRS(lcoeff(q2, v), rs, R) then
    #     qMonic := normal(q2/lcoeff(q2, v));
    #     if denom(qMonic) = 1 then
    #         q2 := qMonic;
    #     end if;
    # end if;
    
    return [q1, q2];
    
end proc;