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