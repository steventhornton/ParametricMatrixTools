# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isEqualOverCS.mpl                                                       #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Determine if two polynomials are equal at all points in a constructible #
# set.                                                                    #
#                                                                         #
# INPUT                                                                   #
#   p1 ... Polynomial                                                     #
#   p2 ... Polynomial                                                     #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   True if p1 = p2 at all points in cs. False otherwise.                 #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([x, a]);                                        #
#   > p1 := 2*x^2 - 12:                                                   #
#   > p2 := (a-1)*x^2 - 4*a:                                              #
#   > cs := GeneralConstruct([a-3], [], R):                               #
#   > isEqualOverCS(p1, p2, cs, R);                                       #
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
isEqualOverCS := proc(p1::polynom, p2::polynom, cs::TRDcs, R::TRDring, $) :: truefalse;
    
    # p1 must be a polynomial in R
    ASSERT(RC:-TRDis_poly(p1, R), "p1 must be a polynomial in R");
    
    # p2 must be a polynomial in R
    ASSERT(RC:-TRDis_poly(p2, R), "p2 must be a polynomial in R");
    
    return isZeroOverCS(p1 - p2, cs, R);
    
end proc;