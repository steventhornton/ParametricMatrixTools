# ======================================================================= #
# ======================================================================= #
#                                                                         #
# JordanForm.mpl                                                          #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 9/2017                                                 #
#                                                                         #
# Compute the Jordan canonical form of a square matrix. If the matrix     #
# contains parameters, the generic result is returned.                    #
#                                                                         #
# INPUT                                                                   #
#   A ... Square matrix with parameters                                   #
#                                                                         #
# OUTPUT                                                                  #
#   The Jordan form of the input matrix.                                  #
#                                                                         #
# EXAMPLE                                                                 #
#   > p1 := (x + 1)^2 * (x^2 + x + 1) * (x + a):                          #
#   > p2 := (x^2 + x + 1) * (x + a):                                      #
#   > C1 := CompanionMatrix(p1, x):                                       #
#   > C2 := CompanionMatrix(p2, x):                                       #
#   > A := DiagonalMatrix([C1, C2]):                                      #
#   > J := JordanForm(A):                                                 #
#   > convert(J, radical);                                                #
#         [-1/2+((1/2)*I)*sqrt(3), 0, 0, 0, 0, 0, 0, 0]                   #
#         [0, -1/2-((1/2)*I)*sqrt(3), 0, 0, 0, 0, 0, 0]                   #
#         [0, 0, -1, 1, 0, 0, 0, 0]                                       #
#         [0, 0, 0, -1, 0, 0, 0, 0]                                       #
#         [0, 0, 0, 0, -a, 0, 0, 0]                                       #
#         [0, 0, 0, 0, 0, -1/2+((1/2)*I)*sqrt(3), 0, 0]                   #
#         [0, 0, 0, 0, 0, 0, -1/2-((1/2)*I)*sqrt(3), 0]                   #
#         [0, 0, 0, 0, 0, 0, 0, -a]                                       #
# ======================================================================= #
# ======================================================================= #
JordanForm := module()

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
#   Same as JordanForm                                                    #
# ----------------------------------------------------------------------- #
init := proc(A::Matrix(square), $)

    # Check the types for the entries of A
    
    return implementation(A);

end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Compute the JCF of A.                                                   #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as JordanFormWith                                                #
# ----------------------------------------------------------------------- #
implementation := proc(A::Matrix(square), $)

    local F :: Matrix,
          J :: Matrix,
          i :: posint,
          blocks,
          charPoly :: polynom,
          charPolyList :: list(polynom),
          q :: polynom,
          m :: posint,
          sqrFreeCharPoly, 
          term :: [polynom, posint];

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

    return J;

end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# EXTERNAL FILES
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveJordanForm/getCharPolys.mpl>

end module;