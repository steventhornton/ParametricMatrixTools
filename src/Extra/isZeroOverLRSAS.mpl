# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isZeroOverLRSAS.mpl                                                     #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Sept. 29/2017                                               #
#                                                                         #
# Determine if a polynomial vanishes at all points in a list of regular   #
# semi-algebraic systems.                                                 #
#                                                                         #
# INPUT                                                                   #
#   p ....... Polynomial                                                  #
#   lrsas ... List of regular semi-algebraic systems                      #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#   - True if p vanishes at all points in the input list of regular       #
#     semi-algebraic systems                                              #
#   - False if there exists a point in the input list of regular          #
#     semi-algebraic systems where p does not vanish                      #
#                                                                         #
# EXAMPLE                                                                 #
#                                                                         #
# ======================================================================= #
# ======================================================================= #
isZeroOverLRSAS := proc(in_p::depends(polyInRing(R)), lrsas::TRDlrsas, R::TRDring, $) :: truefalse;
    
    local p :: polynom, 
          rsas :: TRDrsas;
    
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
    
    # p must be zero over all representing regular semi-algebraic systems in lrsas
    for rsas in lrsas do
        if not isZeroOverRSAS(p, rsas, R) then
            return false;
        end if;
    end do;
    
    return true;

end proc;