# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isZeroOverRS.mpl                                                        #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 24/2016                                                #
#                                                                         #
# Determine if a polynomial vanishes everywhere in a regular system. That #
# is, for a polynomial p, return true if p is in the radical of the       #
# saturated ideal of T where T is a regular chain. The input regular      #
# system represents all points in W(T)\V(h) where h is a polynomial       #
# representing the inequation constraints of the system.                  #
#                                                                         #
# INPUT                                                                   #
#   p .... Polynomial in R                                                #
#   rs ... Regular system of R                                            #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   - True if p vanishes at all points in the input regular system        #
#   - False if a point exists in the input regular system where p does    #
#     not vanish                                                          #
#                                                                         #
# EXAMPLE                                                                 #
#                                                                         #
# ======================================================================= #
# ======================================================================= #
isZeroOverRS := proc(in_p::depends(polyInRing(R)), rs::TRDrs, R::TRDring, $) :: truefalse;
    
    local p :: polynom;
    
    p := RC:-TRDmodularize_coefficients(in_p, R);
    
    # Trivial case where p is a constant
    if RC:-TRDis_constant(p, R) then
        return evalb(p = 0);
    end if;
    
    # Determine if p is in the radical of the saturated ideal of T
    return RC:-TRDis_in_radical(p, RC:-TRDregular_chain_rs(rs, R), R);
    
end proc;