# ======================================================================= #
# ======================================================================= #
#                                                                         #
# TRDis_partiton_cs.mpl                                                   #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 23/2016                                                #
#                                                                         #
# Determine a list of constructible sets forms a partition of another     #
# constructible set.                                                      #
#                                                                         #
# INPUT                                                                   #
#   lcs ... List of constructible sets                                    #
#   cs .... Constructible set                                             #
#   R ..... Polynomial ring                                               #
#                                                                         #
# OUTPUT                                                                  #
#   True if lcs forms a partition of cs, false otherwise.                 #
#   - lcs must be a pairwise disjoint list of constructible sets whose    #
#     union is equal to cs.                                               #
# ======================================================================= #
# ======================================================================= #
TRDis_partition_cs := proc(lcs::TRDlcs, cs::TRDcs, R::TRDring, $) :: truefalse;

    local lcs_nonempty :: TRDlcs,
          csU :: TRDcs;
          
    
    # Remove any empty constructible sets
    lcs_nonempty := remove(c -> RC:-TRDis_empty_constructible_set(c, R), lcs);
    
    # Check that the union of all constructible sets in lcs is equal to cs
    csU := ListUnion(lcs_nonempty, R);
    if not TRDequal_cs(cs, csU, R) then
        return false;
    end if;
    
    return TRDis_disjoint_lcs(lcs_nonempty, R);
    
end proc;