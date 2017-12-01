# ======================================================================= #
# ======================================================================= #
#                                                                         #
# cleanCS.mpl                                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 23/2016                                                #
#                                                                         #
# Given a list with elements of the form                                  #
#   [g, cs]                                                               #
# clean the result by combining cases with the same gcd.                  #
#                                                                         #
# INPUT                                                                   #
#   result ... A list with elements of the form                           #
#                  [g, cs]                                                #
#              where g is a polynomial and rs is a regular system.        #
#   R ........ Polynomial ring                                            #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [g, cs]                                                           #
#   where any cases with the same gcd in the input list have been         #
#   combined.                                                             #
# ======================================================================= #
# ======================================================================= #
cleanCS := proc(result::list([polynom, TRDcs]), R::TRDring, $)

    local csTasks :: Stack,
          item :: [polynom, TRDcs],
          out,
          g_i :: polynom,
          cs_i :: TRDcs,
          noMatch :: truefalse,
          j :: posint,
          g_j :: polynom,
          cs_j :: TRDcs,
          equalAll :: truefalse,
          cs :: TRDcs;

    csTasks := SimpleStack();

    # Add everything from result to the stack
    for item in result do
        csTasks:-push(item);
    end do;

    # Output list
    out := [];

    while not csTasks:-empty() do
        g_i, cs_i := op(csTasks:-pop());

        noMatch := true;

        for j to nops(out) while noMatch do
            g_j, cs_j := op(out[j]);

            equalAll := true;
            if not isZeroOverCS(g_i - g_i, cs_j, R) or not isZeroOverCS(g_i - g_j, cs_i, R) then
                equalAll := false;
            end if;

            # Add to current item in out
            if equalAll then
                noMatch := false;
                cs := RC_CST:-Union(cs_i, cs_j, R);
                out[j] := [g_j, cs];
            end if;

        end do;

        # If no match is found, add to the list
        if noMatch then
            out := [op(out), [g_i, cs_i]];
        end if;

    end do;

    return out;

end proc;