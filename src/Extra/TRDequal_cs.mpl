# ======================================================================= #
# ======================================================================= #
#                                                                         #
# TRDequal_cs.mpl                                                         #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 23/2016                                                #
#                                                                         #
# Determine if two constructible sets are equal.                          #
#                                                                         #
# INPUT                                                                   #
#   cs1 ... Constructible set                                             #
#   cs2 ... Constructible set                                             #
#   R ..... Polynomial ring                                               #
#                                                                         #
# OUTPUT                                                                  #
#   True if cs1 equals cs2, false otherwise.                              #
# ======================================================================= #
# ======================================================================= #
TRDequal_cs := proc(cs1::TRDcs, cs2::TRDcs, R::TRDring, $) :: truefalse;
    
    return RC:-TRDconstructible_set_is_contained(cs1, cs2, R) and RC:-TRDconstructible_set_is_contained(cs2, cs1, R);
    
end proc;