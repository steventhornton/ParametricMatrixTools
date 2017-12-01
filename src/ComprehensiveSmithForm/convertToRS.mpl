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
#   [S, cs]                                                               #
# expand the list by extracting all regular systems from each             #
# constructible set.                                                      #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [S, cs]                                                #
#              where S is a matrix and cs is a constructible set.         #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [S, rs]                                                           #
# ======================================================================= #
# ======================================================================= #
convertToRS := proc(result, R::TRDring, $)

    local output:: {[], list([Matrix, TRDrs])},
          pair :: [Matrix, TRDcs],
          S :: Matrix,
          cs :: TRDcs,
          lrs :: TRDlrs,
          rs :: TRDrs;
    
    output := [];
    
    for pair in result do
        S, cs := op(pair);
        lrs := RC_CST:-RepresentingRegularSystems(cs, R);
        for rs in lrs do
          output := [op(output), [S, rs]];
        end do;
    end do;
    
    return output;

end proc;