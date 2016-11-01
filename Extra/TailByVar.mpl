# ======================================================================= #
# ======================================================================= #
#                                                                         #
# TailByVar.mpl                                                           #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Returns the tail of a polynomial w.r.t a given variable.                #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   TailByVar(p, v, R)                                                    #
#                                                                         #
# INPUT                                                                   #
#   p ... Polynomial                                                      #
#   v ... Variable                                                        #
#   R ... Polynomial ring                                                 #
#                                                                         #
# OUTPUT                                                                  #
#   p - the term containing the largest degree in v.                      #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([x, a, b]):                                     #
#   > p := x^2 + 2x + 1:                                                  #
#   > TailByVar(p, x, R);                                                 #
#         2x + 1                                                          #
#   > p := 10:                                                            #
#   > TailByVar(p, x, R);                                                 #
#         0                                                               #
#   > p := 10a^2 + 13b^3 + a*b - 4:                                       #
#   > TailByVar(p, x, R);                                                 #
#         0                                                               #
#   > p := (a^2 + b)x^3 + (a*b + a^2 - b^2 -1)x - 20:                     #
#   > TailByVar(p, x, R);                                                 #
#         (a*b + a^2 - b^2 -1)x - 20                                      #
#   > p := 0:                                                             #
#   > TailByVar(p, x, R);                                                 #
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
TailByVar := proc(p::polynom, v::name, R::TRDring, $) :: polynom;

    local d :: nonnegint, 
          l::{list(polynom), polynom};

    # p must be a polynomial in R
    ASSERT(RC:-TRDis_poly(p, R), "p must be a polynomial in R");

    if RC:-TRDis_constant(p, R) then
        return 0;
    end if;

    d := degree(p,v);

    if d = 0 then
        return 0;
    end if;

    l := PT:-CoefficientList(p, v);
    l := l[1..-2];
    l := PT:-FromCoefficientList(l, v);

    return l;

end proc;