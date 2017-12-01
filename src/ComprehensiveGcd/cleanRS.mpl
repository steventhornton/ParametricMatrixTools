# ======================================================================= #
# ======================================================================= #
#                                                                         #
# cleanRS.mpl                                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 24/2017                                                #
#                                                                         #
# Given a list with elements of the form                                  #
#   [g, rs]                                                               #
# clean the polynomial g by:                                              #
#    - removing its content                                               #
#    - reducing modulo the regular chain of rs                            #
#    - Ensure the sign(lcoeff(g)) is positive                             #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [g, rs]                                                #
#              where g is a polynomial and rs is a regular system.        #
#   v ........ Variable                                                   #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [g, rs]                                                           #
#   with the same number of elements as the input list.                   #
# ======================================================================= #
# ======================================================================= #
cleanRS := proc(result::list([polynom, TRDrs]), v::name, R::TRDring, $)

    local out :: list([polynom, TRDrs]),
          item :: [polynom, TRDrs],
          g :: polynom,
          rs :: TRDrs,
          rc :: TRDrc;
    
    out := [];
    
    for item in result do
        g, rs := op(item);
        rc := RC_CST:-RepresentingChain(rs, R);
        
        g := primpart(RC:-SparsePseudoRemainder(g, rc, R), v);
        g := sign(lcoeff(g,v))*g;
        
        out := [op(out), [g, rs]];
    end do;
    
    return out;
    
end proc;