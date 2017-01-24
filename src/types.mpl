# ======================================================================= #
# ======================================================================= #
#                                                                         #
# types.mpl                                                               #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 15/2016                                                #
#                                                                         #
# Definitions of types used for type checking in functions. The loadTypes #
# function is called when the ParametricMatrixTools module is loaded.     #
#                                                                         #
# TYPES                                                                   #
#   polyInRing                                                            #
#   TRDring                                                               #
#   TRDrc                                                                 #
#   TRDcs                                                                 #
#   TRDrs                                                                 #
#   TRDsrc                                                                #
#   TRDlrs                                                                #
#   TRDlrc                                                                #
#   TRDlcs                                                                #
# ======================================================================= #
# ======================================================================= #
loadTypes := proc()
    
    userinfo(2, 'ParametricMatrixTools', "Adding types.");
    
    # ------------------------------------------------------------------- #
    # polyInRing                                                          #
    #                                                                     #
    # Check if a polynomial belongs to the given ring.                    #
    # ------------------------------------------------------------------- #
    TypeTools:-AddType('polyInRing', proc(p, R) 
        RC:-TRDis_poly(p, R);
    end proc);
    
    
    # ------------------------------------------------------------------- #
    # TRDring                                                             #
    #                                                                     #
    # Check if R is a valid polynomial ring.                              #
    # ------------------------------------------------------------------- #
    TypeTools:-AddType('TRDring', proc(R)
        RC:-TRDis_polynomial_ring(R);
    end proc);
    
    
    # ------------------------------------------------------------------- #
    # TRDrc                                                               #
    #                                                                     #
    # Check if rc is a valid regular chain.                               #
    # ------------------------------------------------------------------- #
    TypeTools:-AddType('TRDrc', proc(rc)
        RC:-TRDis_regular_chain(rc);
    end proc);
    
    
    # ------------------------------------------------------------------- #
    # TRDcs                                                               #
    #                                                                     #
    # Check if cs is a valid constructible set.                           #
    # ------------------------------------------------------------------- #
    TypeTools:-AddType('TRDcs', proc(cs)
        RC:-TRDis_constructible_set(cs);
    end proc);
    
    
    # ------------------------------------------------------------------- #
    # TRDrs                                                               #
    #                                                                     #
    # Check if rs is a valid regular system.                              #
    # ------------------------------------------------------------------- #
    TypeTools:-AddType('TRDrs', proc(rs)
        RC:-TRDis_regular_system(rs);
    end proc);
    
    
    # ------------------------------------------------------------------- #
    # TRDsrc                                                              #
    #                                                                     #
    # Check if src is a valid subresultant chain.                         #
    # ------------------------------------------------------------------- #
    TypeTools:-AddType('TRDsrc', proc(src)
        RC:-TRDis_subresultant_chain(src);
    end proc);
    
    
    # ------------------------------------------------------------------- #
    # TRDlrs                                                              #
    #                                                                     #
    # Check if lrs is a valid list of regular systems.                    #
    # ------------------------------------------------------------------- #
    TypeTools:-AddType('TRDlrs', proc(lrs)
        type(lrs, list('TRDrs'));
    end proc);
    
    
    # ------------------------------------------------------------------- #
    # TRDlcs                                                              #
    #                                                                     #
    # Check if lcs is a valid list of constructible sets.                 #
    # ------------------------------------------------------------------- #
    TypeTools:-AddType('TRDlcs', proc(lcs)
        type(lcs, list('TRDcs'));
    end proc);
    
    
    # ------------------------------------------------------------------- #
    # TRDlrs                                                              #
    #                                                                     #
    # Check if lrs is a valid list of regular chains.                     #
    # ------------------------------------------------------------------- #
    TypeTools:-AddType('TRDlrc', proc(lrc)
        type(lrc, list('TRDrc'));
    end proc);

end proc: