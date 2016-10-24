# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ParametricSmithForm_Minors.mpl                                          #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 24/2016                                                #
#                                                                         #
# Computes the Smith Normal Form of a matrix where the entries are        #
# univariate polynomials with multivariate coefficients in the parameters #
# by computing gcd's of minors of the input matrix and dividing           #
# successive gcd's. Computation is done modulo a regular system or        #
# constructible set.                                                      #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   ParametricSmithForm_Minors(A, v, cs, R)                               #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   v .... Variable                                                       #
#   cs ... Constructible set                                              #
#   R .... Polynomial Ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [S_i, cs_i]                                                       #
#   where S_i is the Smith Normal Form of A for all parameter values that #
#   satisfy the equations and inequations of cs_i. Together, all the      #
#   constructible sets cs_i (for all values of i) represent partition of  #
#   the input constructible set.                                          #
#                                                                         #
# EXAMPLE                                                                 #
#                                                                         #
# REFERENCES                                                              #
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
ParametricSmithForm_Minors := module()

    export ModuleApply;

    local
        implementation_cs,
        computeDeterminantDivisors,
        determinantDivisor,
        getMinor,
        computeInvariantFactors,
        divideEntries,
        makeDeterminantDivisorsDisjoint,
        convertToMonic;
    
    ModuleApply := proc(A::Matrix(square), v::name, cs::TRDcs, R::TRDring, $) :: list([Matrix, TRDcs]);
        return implementation_cs(A, v, cs, R);
    end proc:

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ----------------------------------------------------------------------- #
# implementation_cs                                                       #
#                                                                         #
# Compute the parametric Smith Normal Form of a matrix by computing gcd's #
# of minors of the input matrix and dividing successive gcd's.            #
#                                                                         #
# INPUT/OUTPUT                                                            #
#    Same as ParametricSmithForm_Minors                                   #
# ----------------------------------------------------------------------- #
implementation_cs := proc(A::Matrix(square), v::name, cs::TRDcs, R::TRDring, $) :: list([Matrix, TRDcs]);
    
    local detDivisors :: 'Array'('datatype'=['Array'('datatype'=polynom), TRDcs]),
          invariantFactors :: 'Array'('datatype'=['Array'('datatype'=ratpoly), TRDcs]);
    
    # An array where each entry is of the form:
    #   [detDiv, cs]
    # where detDiv is an array (index starting at 0) containing n+1 entries
    # where detDiv[i] is the ith determinant divisor of A and is continuous
    # for all parameter values satisfying cs. The collection of all 
    # constructible sets in detDivisors forms a partition of the input
    # constructible set.
    detDivisors := computeDeterminantDivisors(A, v, cs, R);
    
    # Compute the entries of the SNF by dividing successive determinant 
    # minors.
    invariantFactors := computeInvariantFactors(detDivisors, R);
    
    # Construct SNF for each entry of invariantFactors
    return convert(map(input -> [LA:-DiagonalMatrix(map(collect, input[1], v)), 
                                 RC:-TRDrename_constructible_set(input[2])], 
                                 invariantFactors), list);
    
end proc:


# ----------------------------------------------------------------------- #
# computeDeterminantDivisors                                              #
#                                                                         #
# Compute all determinant divisors of a matrix.                           #
#                                                                         #
# INPUT                                                                   #
#    A .... Matrix                                                        #
#    v .... Variable                                                      #
#    cs ... Constructible set                                             #
#    R .... Polynomial Ring                                               #
#                                                                         #
# OUTPUT                                                                  #
#   detDivisors: An array containing n+1 entries (indexed from 0 to n)    #
#   where detDivisors[i] is a list with elements of the form:             #
#       [p, cs]                                                           #
#   where p is the gcd of all i x i minors of A over cs. 0 <= i <= n      #
# ----------------------------------------------------------------------- #
computeDeterminantDivisors := proc(A::Matrix(square), v::name, cs::TRDcs, R::TRDring, $) :: 'Array'('datatype'=['Array'('datatype'=ratpoly), TRDcs]);

    local detDivisors :: 'Array'('datatype'=list([polynom, 'TRDcs'])),
          n :: posint,
          i :: posint;

    n := LA:-RowDimension(A);
    
    detDivisors := Array(0..n, 'datatype'=list([polynom, 'TRDcs']), 'fill'=[[0, RC:-TRDempty_constructible_set()]]);

    # Initialize the first and last elements
    detDivisors[0] := [[1, cs]];
    detDivisors[n] := [[LA:-Determinant(A), cs]];

    for i to n-1 do
        detDivisors[i] := determinantDivisor(A, i, v, cs, R);
    end do;
    
    # Make detDivisors disjoint (Splitting)
    return makeDeterminantDivisorsDisjoint(detDivisors, R);
    
end proc;


# ----------------------------------------------------------------------- #
# determinantDivisor                                                      #
#                                                                         #
# Compute the i-th determinant divisor of a matrix.                       #
#                                                                         #
# INPUT                                                                   #
#    A .... Matrix                                                        #
#    i .... Positive integer                                              #
#    v .... Variable                                                      #
#    cs ... Constructible set                                             #
#    R .... Polynomial ring                                               #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form:                                     #
#       [p, cs_out]                                                       #
#   where p is the ith determinant divisor of A over cs_out.              #
# ----------------------------------------------------------------------- #
determinantDivisor := proc(A::Matrix(square), i::posint, v::name, cs::TRDcs, R::TRDring, $) :: list([polynom, TRDcs]);

    local n :: posint,
          j :: posint,
          k :: posint,
          r :: listlist(posint),
          minors :: 'Array'('datatype'=polynom);

    # Compute all ixi minors
    n := LA:-RowDimension(A);

    ASSERT(i <= n, "Minor matrix size must be less than matrix size");

    r := combinat:-choose(n, n-i);

    minors := Array(1..(nops(r)^2), 'datatype'=polynom);

    # The most expensive part of algorithm!
    for j to nops(r) do
        for k to nops(r) do
            minors[j + nops(r)*(k-1)] := getMinor(A, r[j], r[k]);
        end do;
    end do;
    
    # Compute the gcd of all the minors
    return ListParametricGcd(convert(minors, 'list'), v, cs, R, 'outputType'='CS', 'lazy'=false);

end proc;


# ----------------------------------------------------------------------- #
# getMinor                                                                #
#                                                                         #
# Get the minor of A by removing rows in rL and columns in cL.            #
#                                                                         #
# INPUT                                                                   #
#    A .... Matrix                                                        #
#    rL ... List of positive integers (rows)                              #
#    cL ... List of positive integers (columns)                           #
#                                                                         #
# OUTPUT                                                                  #
#    The determinant of the minor of A by removing rows in rL and columns #
#    in cL.                                                               #
# ----------------------------------------------------------------------- #
getMinor := proc(A::~Matrix(square), rL::list(posint), cL::list(posint), $) :: polynom;

    local B :: Matrix;

    # Ensure rSet and cSet are of the same size
    ASSERT(nops(rL)=nops(cL), "rL and cL must have the same number of elements");

    # Delete the rows and columns of A
    B := LA:-DeleteRow(A, rL);
    B := LA:-DeleteColumn(B, cL);

    return LA:-Determinant(B);

end proc;


# ----------------------------------------------------------------------- #
# computeEntries                                                          #
#                                                                         #
# Compute the entries of the SNF by dividing successive determinant       #
# minors.                                                                 #
#                                                                         #
# INPUT                                                                   #
#   detDivisors ... An Array containing n+1 entries where detDivisors[i]  #
#                   is a list with elements of the form [p, cs] where p   #
#                   is the gcd of all i x i minors of A over cs.          #
#                   0 <= i <= n                                           #
#   v ............. Variable                                              #
#   R ............. Polynomial Ring                                       #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form:                                     #
#        [polyArray, cs]                                                  #
#   where polyArray is an array of rational polynomials and cs is a       #
#   constructible set. The rational polynomials in polyArray are the      #
#   invariant factors of the original matrix.                             #
# ----------------------------------------------------------------------- #
computeInvariantFactors := proc(detDivisors::Array([Array(polynom), TRDcs]), R::TRDring, $) :: 'Array'('datatype'=['Array'('datatype'=ratpoly), TRDcs]);
    
    # Divide successive entries in detDivisors
    return map(input -> [divideEntries(input[1], input[2], R), 
                         input[2]], 
                         detDivisors);
    
end proc:


# ----------------------------------------------------------------------- #
# divideEntries                                                           #
#                                                                         #
# Divide the successive polynomials in an Array by each other.            #
#                                                                         #
# INPUT                                                                   #
#   polyArray ... Array of polynomials (indexed starting at 0)            #
#   cs .......... Constructible set                                       #
#   R ........... Polynomial ring                                         #
#                                                                         #
# OUTPUT                                                                  #
#   A list of polynomials with 1 less element that the input polyArray    #
#   such that                                                             #
#       outPolyList[1] = polyArray[1]/polyArray[0]                        #
#       outPolyList[2] = polyArray[2]/polyArray[1]                        #
#       ...                                                               #
#                                                                         #
# TO DO                                                                   #
#   pseudo-division?                                                      #
# ----------------------------------------------------------------------- #
divideEntries := proc(polyArray::Array(ratpoly), cs::TRDcs, R::TRDring, $) :: 'Array'('datatype'=ratpoly);

    local n :: posint,
          i :: posint,
          result :: 'Array'('datatype'=ratpoly);

    n := ArrayNumElems(polyArray)-1;

    result := Array(1..n, 'datatype'=ratpoly, 'fill'=0);

    for i from 1 to n do

        # If a zero is found, return
        if isZeroOverCS(polyArray[i-1], cs, R) then
            return result;
        else
            result[i] := convertToMonic(normal(polyArray[i]/polyArray[i-1]));
        end if;

    end do;

    return result;

end proc:


# ----------------------------------------------------------------------- #
# makeDeterminantDivisorsDisjoint                                         #
#                                                                         #
# Split up the Array of the i x i minors into a list where each entry     #
# contains an Array simiar to the input array, along with a constructible #
# set where each constructible set in the output array is disjoint.       #
#                                                                         #
# INPUT                                                                   #
#   dList ... An Array containing n + 1 entries where dList[i] is a list  #
#             with elements of the form [p, cs] where p is the gcd of all #
#             i x i minors of A over cs. 0 <= i <= n                      #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form:                                     #
#       [polyArray, cs]                                                   #
#   where polyArray is an Array of polynomials such that polyArray[i] is  #
#   the gcd of all i x i minors of the input matrix. cs is a              #
#   constructible set such that the polynomials in polyArray are the      #
#   minors for all parameter values that satisfy the equations and        #
#   inequations of cs.                                                    #
# ----------------------------------------------------------------------- #
makeDeterminantDivisorsDisjoint := proc(detDivisors::Array(list([polynom, TRDcs])), R::TRDring, $) :: 'Array'('datatype'=['Array'('datatype'=polynom), TRDcs]);

    local n :: posint,
          countNumEntries :: 'Array'('datatype'=list(posint)),
          allCombinations :: list(list(posint)),
          result :: 'Array'('datatype'=['Array'('datatype'=polynom), TRDcs]),
          c :: list(posint),
          polyArray :: 'Array'('datatype'=polynom),
          csArray :: 'Array'('datatype'=TRDcs),
          i :: posint,
          j :: nonnegint,
          csIntersect :: TRDcs;

    # Get the size of the matrix
    n := nops([entries(detDivisors, 'nolist')]) - 1;

    # Make an array where each element is a list counting from 1 to the
    # number of entries in detDivisors[i].
    countNumEntries := map(proc(lp) local j; options operator, arrow; 
                               [seq(j, j = 1..nops(lp))]; 
                           end proc, detDivisors);

    # Make a list with all combinations of the indices
    allCombinations := allComb(convert(countNumEntries, list));

    # The output list
    result := Array(1..nops(allCombinations), 'datatype'=[Array(polynom), 'TRDcs']);

    # Iterate over all combinations and make disjoint
    for i to nops(allCombinations) do
        c := allCombinations[i];

        # Get the Array of polynomials and Array of constructible sets 
        # corresponding to the entries of dList at the indicies in c
        polyArray := Array(0..n, 'datatype'=polynom);
        csArray := Array(0..n, 'datatype'='TRDcs', 'fill'=RC:-TRDempty_constructible_set());
        
        for j from 0 to n do
            polyArray[j], csArray[j] := op(detDivisors[j][c[j+1]]);
        end do;

        # Compute intersection of all constructible sets
        csIntersect := Intersection(convert(csArray, list), R);

        result[i] := [copy(polyArray), csIntersect];

    end do;

    return result;

end proc;


# ----------------------------------------------------------------------- #
# convertToMonic                                                          #
#                                                                         #
# Make a polynomial monic by dividing by the leading coefficient.         #
#                                                                         #
# INPUT                                                                   #
#   p ... Polynomial                                                      #
#   v ... Variable                                                        #
#                                                                         #
# OUTPUT                                                                  #
#    p/lc(p)                                                              #
# ----------------------------------------------------------------------- #
convertToMonic := proc(p::polynom, v::name, $) :: polynom;

    local q :: polynom;

    if nops(indets(p)) = 0 then
        return 1;
    end if;

    q := coeff(p, v, degree(p,v));

    return normal(p/q);

end proc;

end module: