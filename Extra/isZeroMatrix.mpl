# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isZeroMatrix.mpl                                                        #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 24/2016                                                #
#                                                                         #
# Determine if all entries in a matrix are zero.                          #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   isZeroMatrix(A)                                                       #
#   isZeroMatrix(A, rs, R)                                                #
#   isZeroMatrix(A, cs, R)                                                #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   rs ... Regular system                                                 #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   isZeroMatrix(A)        -> true if all entries are indentically zero,  #
#                             false otherwise                             #
#   isZeroMatrix(A, rs, R) -> true if all entries are zero over rs,       #
#                             false otherwise                             #
#   isZeroMatrix(A, cs, R) -> true if all entries are zero over cs,       #
#                             false otherwise                             #
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
isZeroMatrix := overload(
    [
        proc(A::Matrix, $) :: truefalse;
            option overload;
            return LA:-Equal(A, Matrix(LA:-Dimension(A)));
        end,
        
        proc(A::Matrix, rs::TRDrs, R::TRDring, $) :: truefalse;
            
            option overload;
            
            local n :: posint, 
                  m :: posint, 
                  i :: posint, 
                  j :: posint;
            
            n, m := LA:-Dimension(A);
            
            for i to n do
                for j to m do
                    if not isZeroOverRS(A[i,j], rs, R) then
                        return false;
                    end if;
                end do;
            end do;
            
            return true;
            
        end,
        
        proc(A::Matrix, cs::TRDcs, R::TRDring, $) :: truefalse;
            
            option overload;
            
            local n :: posint, 
                  m :: posint, 
                  i :: posint, 
                  j :: posint;
            
            n, m := LA:-Dimension(A);
            
            for i to n do
                for j to m do
                    if not isZeroOverCS(A[i,j], cs, R) then
                        return false;
                    end if;
                end do;
            end do;
            
            return true;
            
        end
    ]
);