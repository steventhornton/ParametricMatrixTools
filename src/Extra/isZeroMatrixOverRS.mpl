# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isZeroMatrixOverRS.mpl                                                  #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 9/2017                                                 #
#                                                                         #
# Determine if all entries of a matrix of polynomials vanish at all       #
# points in the zero set of a regular system.                             #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   rs ... Regular system                                                 #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#     True if all entries in A vanish everywhere over rs.                 #
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
isZeroMatrixOverRS := proc(A::~Matrix, rs::TRDrs, R::TRDring, $) :: truefalse;

    local n::posint, m::posint, i::posint, j::posint;

    n, m := LA:-Dimension(A);

    for i to n do
        for j to m do
            if not isZeroOverRS(A[i,j], rs, R) then
                return false;
            end if;
        end do;
    end do;

    return true;

end proc;
