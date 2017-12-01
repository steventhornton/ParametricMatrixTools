# ======================================================================= #
# ======================================================================= #
#                                                                         #
# pseudo_cofactor.mpl                                                     #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 24/2017                                                #
#                                                                         #
# Computes the polynomial p/g where g is known to be a factor of p over   #
# a regular chain.                                                        #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   g .... A factor (gcd) of p                                            #
#   v .... Variable                                                       #
#   rc ... Regular chain                                                  #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   The rational expression p/g mod rc.                                   #
#                                                                         #
# ASSUMPTIONS                                                             #
#   g is never 0                                                          #
# ======================================================================= #
# ======================================================================= #
pseudo_cofactor := proc(p::depends(polyInRing(R)), g::depends(polyInRing(R)), v::name, rc::TRDrc, R::TRDring, $)
    
    local r :: polynom,
          m :: polynom,
          q :: ratpoly;
    
    ASSERT(not isZeroOverRS(g, RC_CST:-RegularSystem(rc, [], R), R), "g must not be zero");
    
    r := sprem(p, g, v, 'm', 'q');
    
    ASSERT(isZeroOverRS(r, RC_CST:-RegularSystem(rc, R), R), "Remainder must be zero");
    
    q := normal(q/m);
    
    ASSERT(not v in indets(denom(q)), "Denominator must not contain v");
    
    return q;

end proc;