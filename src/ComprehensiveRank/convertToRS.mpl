# ======================================================================= #
# ======================================================================= #
#                                                                         #
# convertToRS.mpl                                                         #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Sept. 21/2017                                               #
#                                                                         #
# Given a list with elements of the form                                  #
#   [r, cs]                                                               #
# expand the list by extracting all regular systems from each             #
# constructible set.                                                      #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [r, cs]                                                #
#              where r is an integer and cs is a constructible set.       #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [r, rs]                                                           #
# ======================================================================= #
# ======================================================================= #
convertToRS := proc(result, R::TRDring, $)

    local output:: {[], list([nonnegint, TRDrs])},
          pair :: [nonnegint, TRDcs],
          r :: nonnegint,
          cs :: TRDcs,
          lrs :: TRDlrs,
          rs :: TRDrs;

    output := [];

    for pair in result do
        r, cs := op(pair);
        lrs := RC_CST:-RepresentingRegularSystems(cs, R);
        for rs in lrs do
          output := [op(output), [r, rs]];
        end do;
    end do;

    return output;

end proc;