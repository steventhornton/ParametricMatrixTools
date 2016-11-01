# ======================================================================= #
# ======================================================================= #
#                                                                         #
# containsParameters.mpl                                                  #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Determine if a parametric univariate polynomial contains parameters.    #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   containsParameters(p, v, R)                                           #
#                                                                         #
# INPUT                                                                   #
#   p ... Polynomial                                                      #
#   v ... Variable                                                        #
#   R ... Polynomial ring                                                 #
#                                                                         #
# OUTPUT                                                                  #
#   - True if p contains indeterminants other than just v                 #
#   - False if p has no indeterminants or v is the only indeterminant     #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([x, a, b]):                                     #
#   > p := x^2 + 2x + 1:                                                  #
#   > containsParameters(p, x, R);                                        #
#         false                                                           #
#   > p := 10:                                                            #
#   > containsParameters(p, x, R);                                        #
#         false                                                           #
#   > p := 10a^2 + 13b^3 + a*b - 4:                                       #
#   > containsParameters(p, x, R);                                        #
#         true                                                            #
#   > p := (a^2 + b)x^3 + (a*b + a^2 - b^2 -1)x - 20:                     #
#   > containsParameters(p, x, R);                                        #
#         true                                                            #
#                                                                         #
# ASSUMPTIONS                                                             #
#   v is the greatest variable of R                                       #
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
containsParameters := proc(p::polynom, v::name, R::TRDring, $) :: truefalse;
    
    # v must be the greatest variable of R
    ASSERT(isGreatestVariable(v, R), "v must be the greatest variable of R");
    
    # p must be a polynomial in R
    ASSERT(RC:-TRDis_poly(p, R), "p must be a polynomial in R");
    
    if nops(indets(p)) > 1 then
        return true;
    elif nops(indets(p)) = 1 then
        if not v in indets(p) then
            return true;
        end if;
    end if;
    
    return false;
    
end proc;