# ======================================================================= #
# ======================================================================= #
#                                                                         #
# cleanRS.mpl                                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 23/2016                                                #
#                                                                         #
# Given a list with elements of the form                                  #
#   [g, rs]                                                               #
# clean the polynomial g by dividing by its leading coefficient w.r.t v.  #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [g, rs]                                                #
#              where g is a polynomial and rs is a regular system.        #
#   v ........ Variable                                                   #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [g, rs]                                                           #
#   with the same number of elements as the input list, and in the same   #
#   order as the input list. g has been divided by it's leading           #
#   coefficient.                                                          #
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
cleanRS := proc(result::list([polynom, TRDrs]), v::name, R::TRDring, $)

    local output,
          pair :: [polynom, TRDrs],
          g :: ratpoly,
          rs :: TRDrs,
          rc :: TRDrc, 
          inv,
          zdiv, 
          gFactors,
          gMonic, 
          rc_inv::TRDrc;

    output := [];

    for pair in result do
        g, rs := op(pair);

        rc := RC_CST:-RepresentingChain(rs, R);

        g := RC:-SparsePseudoRemainder(g, rc, R); 

        g := normal(g/listGcd([coeffs(g, v)]));

        gMonic := normal(g/lcoeff(g, v));

        if denom(gMonic) <> 1 then
            inv, zdiv := op(RC:-Inverse(lcoeff(g, v), rc, R));
            rc_inv := inv[1][3];

            ASSERT(nops(zdiv) = 0, "Must not be any zero divisors");
            ASSERT(nops(inv) = 1, "Inverse must only contain one case.");
            ASSERT(RC_CT:-EqualSaturatedIdeals(rc, rc_inv, R), "Inverse must not lose any cases.");

            g := RC:-SparsePseudoRemainder(g*inv[1][1], rc, R);
            if denom(normal(g/lcoeff(g, v))) = 1 then
                g := normal(g/lcoeff(g, v));
            end if;
        end if;

        g := RC:-SparsePseudoRemainder(g, rc, R); 
        gFactors := factors(g)[2];
        gFactors := remove((x,v) -> not v in indets(x[1]), gFactors, v);
        g := mul(map(x -> x[1]^x[2], gFactors));

        output := [op(output), [g, rs]];

    end do;

    return output;

end proc;