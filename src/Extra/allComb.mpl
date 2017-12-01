# ======================================================================= #
# ======================================================================= #
#                                                                         #
# allComb.mpl                                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Compute all ordered combinations of the elements of a list of lists.    #
#                                                                         #
# INPUT                                                                   #
#    l ... list of lists                                                  #
#                                                                         #
# OUTPUT                                                                  #
#    A list of lists as described above.                                  #
#                                                                         #
# EXAMPLE                                                                 #
#    > l := [[1,2,3], [a,b,c], [x,y]]:                                    #
#    > allComb(l);                                                        #
#          [[1,a,x], [1,a,y], [1,b,x], [1,b,y], [1,c,x], [1,c,y], ...     #
#           [2,a,x], [2,a,y], [2,b,x], [2,b,y], [2,c,x], [2,c,y], ...     #
#           [3,a,x], [3,a,y], [3,b,x], [3,b,y], [3,c,x], [3,c,y]]         #
# ======================================================================= #
# ======================================================================= #
allComb := module()
    
    export ModuleApply;
    
    local addOne;
    
    ModuleApply := proc(l::list(list), $) :: list(list);
        
        local result :: list(list), 
              i :: posint;
        
        result := [[]];
        
        for i to nops(l) do
            result := addOne(result, l[i]);
        end do;
        
        return result;
        
    end proc;
    
    # ------------------------------------------------------------------- #
    # addOne                                                              #
    #                                                                     #
    # Part of allComb                                                     #
    # ------------------------------------------------------------------- #
    addOne := proc(ll::list(list), l::list, $) :: list;

        local result, i, j;

        result := [];

        if nops(ll) = 0 then
            return map(`[]`, op(l));
        end if;

        for j to nops(l) do
            for i to nops(ll) do
                result := [op(result), [op(ll[i]), l[j]]];
            end do;
        end do;

        return result;

    end proc;
    
end module;