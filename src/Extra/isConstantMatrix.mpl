# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isConstantMatrix.mpl                                                    #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 9/2017                                                 #
#                                                                         #
# Determine if all entries in a matrix are constant.                      #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   isConstantMatrix(A, R)                                                #
#   isConstantMatrix(A, v, R)                                             #
#                                                                         #
# INPUT                                                                   #
#   A ... Matrix                                                          #
#   v ... Variable                                                        #
#   R ... Polynomial ring                                                 #
#                                                                         #
# OUTPUT                                                                  #
#   3 Arguments:                                                          #
#       True if the all entries of A either have no indeterminants or     #
#       only have v as an indeterminant, false otherwise.                 #
#   2 Arguments:                                                          #
#       True if the all entries of A have no indeterminants, false        #
#       otherwise.                                                        #
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
isConstantMatrix := module()

    export ModuleApply;

    local
        isConstantMatrix2,
        isConstantMatrix3;

    ModuleApply := proc()
        if nargs = 3 then
            return isConstantMatrix3(args)
        elif nargs = 2 then
            return isConstantMatrix2(args)
        else
            error "Function takes 2 or 3 arguments";
        end if;
    end proc;

    # 3 argument method
    isConstantMatrix3 := proc(A::Matrix, v::name, R::TRDring, $) :: truefalse;

        local n::posint, m::posint, i::posint, j::posint;

        # ERROR CHECKING -----------------------------

        # v must be the greatest variable of R
        if not isGreatestVariable(v, R) then
            error "v must be the greatest variable of R";
        end if;

        n, m := LA:-Dimension(A);

        # All elements of A must be polynomials in R
        for i to n do
            for j to m do
                if not RC:-TRDis_poly(A[i,j], R) then
                    error "Expected a matrix of polynomials in R";
                end if;
            end do;
        end do;

        # CODE ---------------------------------------

        for i to n do
            for j to m do
                if not isConstant(A[i,j], v, R) then
                    return false;
                end if;
            end do;
        end do;

        return true;

    end proc;

    # 2 argument method
    isConstantMatrix2 := proc(A::Matrix, R::TRDring, $) :: truefalse;

        local n::posint, m::posint, i::posint, j::posint;

        # ERROR CHECKING -----------------------------

        n, m := LA:-Dimension(A);

        # All elements of A must be polynomials in R
        for i to n do
            for j to m do
                if not RC:-TRDis_poly(A[i,j], R) then
                    error "Expected a matrix of polynomials in R";
                end if;
            end do;
        end do;

        # CODE ---------------------------------------

        for i to n do
            for j to m do
                if not isConstant(A[i,j], R) then
                    return false;
                end if;
            end do;
        end do;

        return true;

    end proc;

end module;