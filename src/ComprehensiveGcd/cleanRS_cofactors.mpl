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
# ======================================================================= #
# ======================================================================= #
cleanRS_cofactors := proc(result::list([polynom, ratpoly, ratpoly, TRDrs]), v::name, R::TRDring, $)

    local item :: [polynom, ratpoly, ratpoly, TRDrs],
          g :: ratpoly,
          c1 :: ratpoly,
          c2 :: ratpoly,
          rs :: TRDrs,
          rc :: TRDrc,
          co :: ratpoly,
          out :: list([polynom, ratpoly, ratpoly, TRDrs]),
          s :: integer,
          h :: polynom,
          h1, h2;

    out := [];

    for item in result do
        g, c1, c2, rs := op(item);
        
        rc := RC_CST:-RepresentingChain(rs, R);
        
        g := RC:-SparsePseudoRemainder(primpart(g, v, 'co'), rc, R, 'h');
        s := sign(lcoeff(g,v));
        g := s*g;

        c1 := normal(c1*co/h);
        c2 := normal(c2*co/h);
        
        c1 := RC:-SparsePseudoRemainder(numer(c1), rc, R, 'h1')/RC:-SparsePseudoRemainder(denom(c1), rc, R, 'h2');
        c1 := normal(c1*h2/h1);
        c2 := RC:-SparsePseudoRemainder(numer(c2), rc, R, 'h1')/RC:-SparsePseudoRemainder(denom(c2), rc, R, 'h2');
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