# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isConstantMatrix.mpl                                                    #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Determine if all entries in a matrix are constant.                      #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   isConstantMatrix(A)                                                   #
#                                                                         #
# INPUT                                                                   #
#   A ... Matrix                                                          #
#                                                                         #
# OUTPUT                                                                  #
#   true if there are no indeterminants in A false otherwise              #
#                                                                         #
# EXAMPLE                                                                 #
#   > A := Matrix(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9]):                     #
#   > isConstantMatrix(A):                                                #
#         true                                                            #
#   > A := Matrix(3, 3, [1, 2, 3, 4, 5, 6, x, y, z]):                     #
#   > isConstantMatrix(A):                                                #
#         false                                                           #
#   > isConstantMatrix(A):                                                #
#         false                                                           #
#   > A := Matrix(3, 4, [x, 2, 3, 4, x, 6, 7, 8, x, 10, 11, 12]):         #
#   > isConstantMatrix(A):                                                #
#         false                                                           #
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
isConstantMatrix := proc(A::Matrix, $) :: truefalse;
    return evalb(nops(indets(A)) = 0);
end proc;