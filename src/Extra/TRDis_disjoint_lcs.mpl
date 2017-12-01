# ======================================================================= #
# ======================================================================= #
#                                                                         #
# TRDis_disjoint_lcs.mpl                                                  #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 23/2016                                                #
#                                                                         #
# Determine if the constructible sets if a list of constructible sets are #
# pairwise disjoint.                                                      #
#                                                                         #
# INPUT                                                                   #
#   lcs ... List of constructible sets                                    #
#   R ..... Polynomial ring                                               #
#                                                                         #
# OUTPUT                                                                  #
#   True if the constructible sets in lcs are pairwise disjoint, false    #
#   otherwise.                                                            #
# ======================================================================= #
# ======================================================================= #
TRDis_disjoint_lcs := proc(lcs::TRDlcs, R::TRDring, $) :: truefalse;

    local lcs_nonempty :: TRDlcs,
          csI :: TRDcs,
          i :: posint,
          j :: posint;
    
    # Remove any empty constructible sets
    lcs_nonempty := remove(c -> RC:-TRDis_empty_constructible_set(c, R), lcs);

    # Check that each pairwise element of lcs_nonempty is disjoint
    if nops(lcs_nonempty) = 1 then
        return true;
    end if;
    for i to nops(lcs_nonempty)-1 do
        for j from i+1 to nops(lcs_nonempty) do
            csI := RC_CST:-Intersection(lcs_nonempty[i], lcs_nonempty[j], R);
            if not RC:-TRDis_empty_constructible_set(csI, R) then
                return false;
            end if;
        end do;
    end do;

    return true;

end proc;