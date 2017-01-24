# ======================================================================= #
# ======================================================================= #
#                                                                         #
# cleanRS_cofactors.mpl                                                   #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 24/2017                                                #
#                                                                         #
# Given a list with elements of the form                                  #
#   [g, c1, c2, rs]                                                       #
# clean the polynomial g by:                                              #
#    - removing its content                                               #
#    - reducing modulo the regular chain of rs                            #
#    - Ensure the sign(lcoeff(g)) is positive                             #
# And ensure that c1*g = p1, and c2*g = p2.                               #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [g, c1, c2, rs]                                        #
#              where g is a polynomial, c1 and c2 are polynomials in v    #
#              where the coefficients are rational functions in the       #
#              remaining variables (parameters) and rs is a regular       #
#              system.                                                    #
#   v ........ Variable                                                   #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [g, c1, c2, rs]                                                   #
#   with the same number of elements as the input list.                   #
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
cleanRS_cofactors := proc(result::list([polynom, ratpoly, ratpoly, TRDrs]), v::name, R::TRDring, $)

    local item :: [polynom, ratpoly, ratpoly, TRDrs],
          g :: polynom,
          c1 :: ratpoly,
          c2 :: ratpoly,
          rs :: TRDrs,
          rc :: TRDrc,
          co :: ratpoly,
          out :: list([polynom, ratpoly, ratpoly, TRDrs]),
          s :: integer,
          h :: polynom,
          h1 :: polynom,
          h2 :: polynom;

    out := [];

    for item in result do
        g, c1, c2, rs := op(item);
        rc := RC_CST:-RepresentingChain(rs, R);

        g := primpart(RC:-SparsePseudoRemainder(g, rc, R, 'h'), v, 'co');
        s := sign(lcoeff(g,v));
        g := s*g;

        c1 := normal(c1*co/h);
        c2 := normal(c2*co/h);

        c1 := normal(RC:-SparsePseudoRemainder(numer(c1), rc, R, 'h1')/RC:-SparsePseudoRemainder(denom(c1), rc, R, 'h2'));
        c1 := normal(c1*h2/h1);
        c2 := normal(RC:-SparsePseudoRemainder(numer(c2), rc, R, 'h1')/RC:-SparsePseudoRemainder(denom(c2), rc, R, 'h2'));
        c2 := normal(c2*h2/h1);

        c1 := s*c1;
        c2 := s*c2;

        # TO DO:
        #   Add heuristic to clean the result by inverting denom(c1) and denom(c2)
        #   over rc.

        out := [op(out), [g, c1, c2, rs]];

    end do;
    
    return out;

end proc;