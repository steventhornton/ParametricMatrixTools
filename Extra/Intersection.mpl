# ======================================================================= #
# ======================================================================= #
#                                                                         #
# Intersection.mpl                                                        #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 24/2016                                                #
#                                                                         #
# Computes the intersection of all elements in a list or set of           #
# constructible sets.                                                     #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   Intersection(csList, R)                                               #
#   Intersection(csSet, R)                                                #
#   Intersection(csArray, R)                                              #
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
#                                                                         #
# LICENSE                                                                 #
#   This program is free software: you can redistribute it and/or modify  #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation, either version 3 of the License, or     #
#   any later version.                                                    #
#                                                                         #
#   This program is distributed in the hope that it will be useful,       #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
#   You should have received a copy of the GNU General Public License     #
#   along with this program.  If not, see http://www.gnu.org/licenses/.   #
# ======================================================================= #
# ======================================================================= #
Intersection := overload(
    [
        proc(csList::list(TRDcs), R::TRDring, $) :: TRDcs; 
            option overload;
            userinfo(2, 'ParametricMatrixTools', "Calling Intersection with a list of constructible sets.");
            return Intersection(convert(csList, Array), R);
        end,

        proc(csSet::set(TRDcs), R::TRDring, $) :: TRDcs; 
            option overload;
            userinfo(2, 'ParametricMatrixTools', "Calling Intersection with a set of constructible sets.");
            return Intersection(convert(csSet, Array), R);
        end,
        
        
        proc(csArray::Array(TRDcs), R::TRDring, $) :: TRDcs; 
            
            option overload;
            
            local n :: nonnegint,
                  i :: posint,
                  cs :: TRDcs;
            
            userinfo(2, 'ParametricMatrixTools', "Calling Intersection with an array of constructible sets.");
            
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