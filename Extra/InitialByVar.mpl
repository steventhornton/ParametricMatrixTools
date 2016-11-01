# ======================================================================= #
# ======================================================================= #
#                                                                         #
# InitialByVar.mpl                                                        #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Returns the initial of a polynomial w.r.t a given variable.             #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   InitialByVar(p, v, R)                                                 #
#                                                                         #
# INPUT                                                                   #
#   p ... Polynomial                                                      #
#   v ... Variable                                                        #
#   R ... Polynomial ring                                                 #
#                                                                         #
# OUTPUT                                                                  #
#   The coefficient in p of the term of largest degree in v.              #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([x, a, b]):                                     #
#   > p := x^2 + 2x + 1:                                                  #
#   > InitialByVar(p, x, R);                                              #
#         1                                                               #
#   > p := 10:                                                            #
#   > InitialByVar(p, x, R);                                              #
#         10                                                              #
#   > p := 10a^2 + 13b^3 + a*b - 4:                                       #
#   > InitialByVar(p, x, R);                                              #
#         10a^2 + 13b^3 + a*b - 4                                         #
#   > p := (a^2 + b)x^3 + (a*b + a^2 - b^2 -1)x - 20:                     #
#   > InitialByVar(p, x, R);                                              #
#         a^2 + b                                                         #
#   > p := 0:                                                             #
#   > InitialByVar(p, x, R);                                              #
#         0                                                               #
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
InitialByVar := proc(p::polynom, v::name, R::TRDring, $) :: polynom;
    
    # p must be a polynomial in R
    ASSERT(RC:-TRDis_poly(p, R), "p must be a polynomial in R");
    
    return coeff(p, v, degree(p, v));
    
end proc;