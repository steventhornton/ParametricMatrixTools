# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isZeroOverRSAS.mpl                                                      #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Sept. 29/2017                                               #
#                                                                         #
# Determine if a polynomial vanishes everywhere in a regular              #
# semi-algebraic system.                                                  #
#                                                                         #
# INPUT                                                                   #
#   p ...... Polynomial in R                                              #
#   rsas ... Regular system of R                                          #
#   R ...... Polynomial ring                                              #
#                                                                         #
# OUTPUT                                                                  #
#   - True if p vanishes at all points in the input regular               #
#     semi-algebraic system                                               #
#   - False if a point exists in the input regular semi-algebraic system  #
#     where p does not vanish                                             #
#                                                                         #
# EXAMPLE                                                                 #
#                                                                         #
# ======================================================================= #
# ======================================================================= #
isZeroOverRSAS := proc(in_p::depends(polyInRing(R)), rsas::TRDrsas, R::TRDring, $) :: truefalse;

    local p :: polynom;

    p := RC:-TRDmodularize_coefficients(in_p, R);

    # Trivial case where p is a constant
    if RC:-TRDis_constant(p, R) then
        return evalb(p = 0);
    end if;

    # Determine if p is in the radical of the saturated ideal of T
    return RC:-TRDis_in_radical(p, RC_SAST:-RepresentingChain(rsas, R), R);

end proc;