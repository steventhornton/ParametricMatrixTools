# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ListIntersection.mpl                                                    #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Computes the intersection of all elements in a list or set of           #
# constructible sets.                                                     #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   ListIntersection(csList, R)                                           #
#   ListIntersection(csSet, R)                                            #
#   ListIntersection(csArray, R)                                          #
#                                                                         #
# INPUT                                                                   #
#   csList .... List or constructible sets                                #
#   csSet ..... Set or constructible sets                                 #
#   csArray ... Array or constructible sets                               #
#   R ......... Polynomial ring                                           #
#                                                                         #
# OUTPUT                                                                  #
#   A constructible set representing the intersection of all input sets   #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([a,b,c]):                                       #
#   >                                                                     #
#   > cs1 := GeneralConstruct([], [], R):                                 #
#   > cs2 := GeneralConstruct([a-1], [], R):                              #
#   > cs3 := GeneralConstruct([a-2], [b-12], R):                          #
#   > cs4 := GeneralConstruct([c+3, b], [b-2, a], R):                     #
#   >                                                                     #
#   > cs5 := GeneralConstruct([a-1, a-2, c+3, b], [b-12, b-2, a], R):     #
#   >                                                                     #
#   > cs := ListIntersection([cs1, cs2, cs3, cs4], R):                    #
#   >                                                                     #
#   > # If cs5 is in cs and cs is in cs5, then they are equal.            #
#   > IsContained(cs5, cs, R) and IsContained(cs, cs5, R);                #
#         true                                                            #
# ======================================================================= #
# ======================================================================= #
ListIntersection := overload(
    [
        proc(csList::list(TRDcs), R::TRDring, $) :: TRDcs; 
            option overload;
            userinfo(2, 'ParametricMatrixTools', "Calling ListIntersection with a list of constructible sets.");
            return ListIntersection(convert(csList, Array), R);
        end,

        proc(csSet::set(TRDcs), R::TRDring, $) :: TRDcs; 
            option overload;
            userinfo(2, 'ParametricMatrixTools', "Calling ListIntersection with a set of constructible sets.");
            return ListIntersection(convert(csSet, Array), R);
        end,
        
        
        proc(csArray::Array, R::TRDring, $) :: TRDcs; 
            
            option overload;
            
            local n :: nonnegint,
                  i :: posint,
                  cs :: TRDcs;
            
            userinfo(2, 'ParametricMatrixTools', "Calling ListIntersection with an array of constructible sets.");
            
            n := ArrayNumElems(csArray);
            
            if n = 1 then
                return csArray[1];
            end if;
            
            # If any of the constructible sets are empty, the intersection
            # of all constructible sets must be empty.
            for cs in csArray do
                if RC_CST:-IsEmpty(cs, R) then
                    return RC:-TRDempty_constructible_set();
                end if;
            end do;
            
            cs := csArray[1];
            
            for i from 2 to n do
                cs := RC_CST:-Intersection(cs, csArray[i], R);
            end do;
            
            return cs;
        end
    ]
);