# ======================================================================= #
# ======================================================================= #
#                                                                         #
# JordanFormWithProviso.mpl                                               #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 21/2016                                                #
#                                                                         #
# Given a matrix A where the entries are multivariate polynomials in      #
# paramters, this method computes the Jordan canonical form of A and      #
# provides conditions on the parameters such that the correct Jordan      #
# form is returned.                                                       #
#                                                                         #
# INPUT                                                                   #
#   A ... Square matrix with parameters                                   #
#                                                                         #
# OUTPUT                                                                  #
#   Two values:                                                           #
#       J, h                                                              #
#   Where J is the Jordan canonical form of A under the constraints       #
#   imposed on the parameters by the inequation h.                        #
#                                                                         #
# EXAMPLE                                                                 #
#   > p1 := (x + 1)^2 * (x^2 + x + 1) * (x + a):                          #
#   > p2 := (x^2 + x + 1) * (x + a):                                      #
#   > C1 := CompanionMatrix(p1, x):                                       #
#   > C2 := CompanionMatrix(p2, x):                                       #
#   > A := DiagonalMatrix([C1, C2]):                                      #
#   > J, h := JordanFormWithProviso(A):                                   #
#   > convert(J, radical);                                                #
#         [-1/2+((1/2)*I)*sqrt(3), 0, 0, 0, 0, 0, 0, 0]                   #
#         [0, -1/2-((1/2)*I)*sqrt(3), 0, 0, 0, 0, 0, 0]                   #
#         [0, 0, -1, 1, 0, 0, 0, 0]                                       #
#         [0, 0, 0, -1, 0, 0, 0, 0]                                       #
#         [0, 0, 0, 0, -a, 0, 0, 0]                                       #
#         [0, 0, 0, 0, 0, -1/2+((1/2)*I)*sqrt(3), 0, 0]                   #
#         [0, 0, 0, 0, 0, 0, -1/2-((1/2)*I)*sqrt(3), 0]                   #
#         [0, 0, 0, 0, 0, 0, 0, -a]                                       #
#   > h;                                                                  #
#         a^3-2*a^2+2*a-1 <> 0                                            #
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
JordanFormWithProviso := module()

    export ModuleApply;

    local init,
          implementation,
          getCharPolys;

    ModuleApply := proc()
        return init(args);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# ----------------------------------------------------------------------- #
# init                                                                    #
#                                                                         #
# Checks the types of the input and calls the implementation if all input #
# values pass checks.                                                     #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as JordanFormWithProviso                                         #
# ----------------------------------------------------------------------- #
init := proc(A::Matrix(square), $)

    # Check if A has parameters
    #if nops(indets(A)) = 0 then
    #    error "Matrix has no parameters, use LinearAlgebra:-JordanForm";
    #end if;

    return implementation(A);

end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Compute the JCF of A.                                                   #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as JordanFormWithProviso                                         #
# ----------------------------------------------------------------------- #
implementation := proc(A::Matrix(square), $)

    local p :: polynom,
          F :: Matrix,
          J :: Matrix,
          i :: posint,
          blocks,
          charPoly :: polynom,
          charPolyList :: list(polynom),
          q :: polynom,
          m :: posint,
          d :: polynom,
          sqrFreeCharPoly, 
          term :: [polynom, posint];

    # Compute the characteristic polynomial of A and compute square-free 
    # factorization.
    p := LA:-CharacteristicPolynomial(A, 'x');
    p := mul(map(y -> y[1], sqrfree(p, 'x')[2]));

    d := discrim(p, 'x');
    d := mul(map(y -> y[1], sqrfree(d)[2]));

    F := LA:-FrobeniusForm(A);

    charPolyList := getCharPolys(F, 'x');

    blocks := [];

    for charPoly in charPolyList do
        sqrFreeCharPoly := sqrfree(factor(charPoly), 'x')[2];
        for term in sqrFreeCharPoly do
            
            q, m := op(term);
            
            if degree(q, 'x') > 1 then
                for i to degree(q, 'x') do
                    blocks := [op(blocks), [RootOf(q, 'x', 'index'=i), m]];
                end do;
            else
                blocks := [op(blocks), [RootOf(q, 'x'), m]];
            end if;
        end do;
    end do;

    J := LA:-JordanBlockMatrix(blocks);

    return J, d <> 0;


end proc;


# ----------------------------------------------------------------------- #
# getCharPolys                                                            #
#                                                                         #
# Extract the characteristic polynomials from a matrix in Frobenius form  #
#                                                                         #
# INPUT                                                                   #
#   F ... Matrix in Frobenius (Rational) normal form                      #
#   v ... The varible to use for the ouput polynomials                    #
#                                                                         #
# OUTPUT                                                                  #
#   A list of the characteristic polynomials from each companion matrix   #
#   in the input matrix.                                                  #
# ----------------------------------------------------------------------- #
getCharPolys := proc(F::Matrix(square), v::name, $) :: list(polynom);

    local subDiag :: list,
          startBlock :: list(posint),
          endBlock :: list(posint),
          i :: posint,
          charPolyList :: list(polynom),
          row1 :: posint, 
          row2 :: posint,
          col :: posint,
          p :: polynom,
          n :: posint;

    n := LA:-RowDimension(F);

    subDiag := convert(LA:-Diagonal(F, -1), list);

    startBlock := [1];
    for i to nops(subDiag) do
        if subDiag[i] = 0 then
            startBlock := [op(startBlock), i+1]
        end if;
    end do;

    endBlock := [op(map(a -> a-1, startBlock[2..-1])), n];

    charPolyList := [];
    for i to nops(startBlock) do
        row1, row2 := startBlock[i], endBlock[i];
        col := endBlock[i];
        p := PT:-FromCoefficientList([op(convert(-F[row1 .. row2, col], list)), 1], v); 
        charPolyList := [op(charPolyList), p];
    end do; 

    return charPolyList;

end proc;

end module;