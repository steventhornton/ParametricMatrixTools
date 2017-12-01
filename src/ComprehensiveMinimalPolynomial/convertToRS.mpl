# ======================================================================= #
# ======================================================================= #
#                                                                         #
# convertToRS.mpl                                                         #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Sept. 30/2017                                               #
#                                                                         #
# Given a list with elements of the form                                  #
#   [pmin, cs]                                                            #
# expand the list by extracting all regular systems from each             #
# constructible set.                                                      #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [pmin, cs]                                             #
#              where pmin is a polynomial and cs is a constructible set.  #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [pmin, rs]                                                        #
# ======================================================================= #
# ======================================================================= #
convertToRS := proc(result, R::TRDring, $)

    local output:: {[], list([polynom, TRDrs])},
          pair :: [polynom, TRDcs],
          pmin :: polynom,
          cs :: TRDcs,
          lrs :: TRDlrs,
          rs :: TRDrs;

    output := [];

    for pair in result do
        pmin, cs := op(pair);
        lrs := RC_CST:-RepresentingRegularSystems(cs, R);
        for rs in lrs do
          output := [op(output), [pmin, rs]];
        end do;
    end do;

    return output;

end proc;