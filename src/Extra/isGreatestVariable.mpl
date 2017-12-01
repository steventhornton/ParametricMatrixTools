# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isGreatestVariable.mpl                                                  #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Determine if a variable is the largest variable in a polynomial ring.   #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   isGreatestVariable(v, R)                                              #
#                                                                         #
# INPUT                                                                   #
#   v ... variable                                                        #
#   R ... Polynomial ring                                                 #
#                                                                         #
# OUTPUT                                                                  #
#   True if v is the largest variable of R, false otherwise.              #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([x, a, b]):                                     #
#   > isGreatestVariable(x, R);                                           #
#         true                                                            #
#   > isGreatestVariable(b, R);                                           #
#         false                                                           #
# ======================================================================= #
# ======================================================================= #
isGreatestVariable := proc(v::name, R::TRDring, $) :: truefalse;
    
    # Ensure v is a variable of R
    ASSERT(evalb(v in R['variables']), "v is not a variable of R");
    
    return evalb(R['variables'][1] = v);
    
end proc;