# ======================================================================= #
# ======================================================================= #
#                                                                         #
# TRDdifference_intersect_cs_p.mpl                                        #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 6/2016                                                 #
#                                                                         #
# Compute the difference and intersection of the zero set of a            #
# constructible set and the zero set of a polynomial.                     #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   TRDdifference_intersect_cs_p(cs, p, R)                                #
#                                                                         #
# INPUT                                                                   #
#   cs ... Constructible set                                              #
#   p .... Polynomial                                                     #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   Two constructible sets cs_nz, cs_z are returned such that p vanishes  #
#   at all points in the zero set of cs_z and p vanishes nowhere in the   #
#   zero set of cs_nz.                                                    #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([a, b]):                                        #
#   > p := a - 1:                                                         #
#   > cs := GeneralConstruct([b + 1], [], R):                             #
#   > cs_nz, cs_z := TRDdifference_intersect_cs_p(cs, p, R):              #
#   > Display(cs_nz, R);                                                  #
#         { b + 1 = 0                                                     #
#         { a - 1 <> 0                                                    #
#   > Display(cs_z, R);                                                   #
#         { a - 1 = 0                                                     #
#         { b - 1 = 0                                                     #
# ======================================================================= #
# ======================================================================= #
TRDdifference_intersect_cs_p := proc(cs::TRDcs, p::depends(polyInRing(R)), R::TRDring, $) :: TRDcs, TRDcs;
    
    local cs_p::TRDcs;
    
    cs_p := RC_CST:-GeneralConstruct([p], [], R);
    
    return RC:-TRDdifference_intersect_cs_cs(cs, cs_p, R);
    
end proc;