# ======================================================================= #
# ======================================================================= #
#                                                                         #
# matrixContainsParameters.mpl                                            #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Determine if a matrix of parametric univariate polynomial contains      #
# parameters.                                                             #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   matrixContainsParameters(A, v, R)                                     #
#                                                                         #
# INPUT                                                                   #
#   A ... Matrix                                                          #
#   v ... Variable                                                        #
#   R ... Polynomial ring                                                 #
#                                                                         #
# OUTPUT                                                                  #
#   - True if an entry exists in A that has an indeterminant other than v #
#   - False otherwise.                                                    #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([x, a, b]):                                     #
#   > A := Matrix(4, 5):                                                  #
#   > matrixContainsParameters(A, x, R);                                  #
#         false                                                           #
#   > A := RandomMatrix(4) + x*RandomMatrix(4):                           #
#   > matrixContainsParameters(A, x, R);                                  #
#         false                                                           #
#   > A := a*RandomMatrix(4) + b*RandomMatrix(4):                         #
#   > matrixContainsParameters(A, x, R);                                  #
#         true                                                            #
#   > A := Matrix(4, 5):                                                  #
#   > A[2, 2] := x:                                                       #
#   > A[3, 2] := a+4*b-2*x:                                               #
#   > matrixContainsParameters(A, x, R);                                  #
#         true                                                            #
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
matrixContainsParameters := proc(A::Matrix, v::name, R::TRDring, $) :: truefalse;

    local i::posint,
          j::posint;

    # v must be the greatest variable of R
    ASSERT(isGreatestVariable(v, R), "v must be the greatest variable of R");
    
    for i to LA:-RowDimension(A) do
        for j to LA:-ColumnDimension(A) do
            if containsParameters(A[i,j], v, R) then
                return true;
            end if;
        end do;
    end do;
    
    return false;
    
end proc;