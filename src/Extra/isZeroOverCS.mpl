# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isZeroOverCS.mpl                                                        #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 24/2016                                                #
#                                                                         #
# Determine if a polynomial vanishes at all points in a constructible     #
# set.                                                                    #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial                                                     #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   - True if p vanishes at all points in the input constructible set     #
#   - False if there exists a point in the input constructible set where  #
#     p does not vanish                                                   #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([a]):                                           #
#   > cs := GeneralConstruct([a-1], [], R):                               #
#   > isZeroOverCS(a+1, cs, R);                                           #
#         false                                                           #
#   > isZeroOverCS(a-1, cs, R);                                           #
#         true                                                            #
# ======================================================================= #
# ======================================================================= #
isZeroOverCS := proc(in_p::depends(polyInRing(R)), cs::TRDcs, R::TRDring, $) :: truefalse;
    
    local p :: polynom, 
          rs :: TRDrs;
    
    # Ensure p is a polynomial in R
    if not RC:-TRDis_poly(in_p, R) then
        error "invalid polynomial: %1", in_p;
    else
        p := RC:-TRDmodularize_coefficients(in_p, R);
    end if;
    
    # Trivial case where p is a constant
    if RC:-TRDis_constant(p, R) then
        return evalb(p = 0);
    end if;
    
    # p must be zero over all representing regular systems in cs
    for rs in RC:-TRDregular_systems(cs, R) do
        if not isZeroOverRS(p, rs, R) then
            return false;
        end if;
    end do;
    
    return true;

end proc;