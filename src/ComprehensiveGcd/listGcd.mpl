# ======================================================================= #
# ======================================================================= #
#                                                                         #
# listGcd.mpl                                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 14/2016                                                #
#                                                                         #
# Compute the gcd of a list of polynomials with one or zero               #
# indeterminants. This is the non-parametric method.                      #
#                                                                         #
# INPUT                                                                   #
#   lp ... List of polynomials                                            #
#                                                                         #
# OUTPUT                                                                  #
#   The gcd of all polynomials in lp.                                     #
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
listGcd := proc(lp::list(polynom), $) :: polynom;

    local g :: polynom,
          i :: posint;

    # Ensure l has at least 2 elements
    if nops(lp) = 1 then
        return lp[1];
    end if;
    
    if nops(lp) = 0 then
        error "lp must contain at least 1 polynomial";
    end if;

    g := gcd(lp[1], lp[2]);

    for i from 3 to nops(lp) do
        g := gcd(lp[i], g);
    end do;

    return g;

end proc;