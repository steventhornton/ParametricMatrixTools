# ======================================================================= #
# ======================================================================= #
#                                                                         #
# getCharPolys.mpl                                                        #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 10/2017                                                #
#                                                                         #
# Extract the characteristic polynomials from a matrix in Frobenius form. #
#                                                                         #
# INPUT                                                                   #
#   F ... Matrix in Frobenius (Rational) normal form                      #
#   v ... The varible to use for the ouput polynomials                    #
#                                                                         #
# OUTPUT                                                                  #
#   A list of the characteristic polynomials from each companion matrix   #
#   in the input matrix.                                                  #
# ======================================================================= #
# ======================================================================= #
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