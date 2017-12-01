# ======================================================================= #
# ======================================================================= #
#                                                                         #
# convertToRS.mpl                                                         #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Sept. 12/2017                                               #
#                                                                         #
# Given a list with elements of the form                                  #
#   [g, cs]                                                               #
# expand the list by extracting all regular systems from each             #
# constructible set.                                                      #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [g, cs]                                                #
#              where g is a polynomial and cs is a constructible set.     #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [g, rs]                                                           #
# ======================================================================= #
# ======================================================================= #
convertToRS := proc(result, R::TRDring, $)

    local output:: {[], list([polynom, TRDrs])},
          pair :: [polynom, TRDcs],
          g :: polynom,
          cs :: TRDcs,
          lrs :: TRDlrs;

    output := [];

    for pair in result do
        g, cs := op(pair);
        lrs := RC_CST:-RepresentingRegularSystems(cs, R);
        output := [op(output), op(zip((x, y) -> [x, y], g, lrs))];
    end do;

    return output;

end proc;