# ======================================================================= #
# ======================================================================= #
#                                                                         #
# convertListWithCSToListWithRS.mpl                                       
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 1/2017                                                 #
#                                                                         #
# Given a list with elements of the form                                  #
#   [..., cs, ...]                                                        #
# expand the list by extracting all regular systems from each             #
# constructible set.                                                      #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [..., cs, ...]                                         #
#              where F is a matrix and cs is a constructible set.         #
#   i ........ Index where the constructible set occurs.                  #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [..., rs, ...]                                                    #
# ======================================================================= #
# ======================================================================= #
convertListWithCSToListWithRS := proc(L::list, i::nonnegint, R::TRDring, $)

    local output :: list(list),
          l1 :: list,
          l2 :: list,
          l :: list,
          lrs :: TRDlrs,
          rs :: TRDrs;

    output := [];
    
    for l in L do
      l1 := l[1..(i-1)];
      l2 := l[(i+1)..nops(l)];
      
      lrs := RC_CST:-RepresentingRegularSystems(l[i], R);
      
      for rs in lrs do
        output := [op(output), [op(l1), rs, op(l2)]];
      end do;
      
    end do;

    return output;

end proc;