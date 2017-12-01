# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isZeroMatrixOverLRSAS.mpl                                               #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Sept. 29/2017                                               #
#                                                                         #
# Determine if all entries of a matrix of polynomials vanish at all       #
# points in the zero set of a list of regular semi-algebraic systems.     #
#                                                                         #
# INPUT                                                                   #
#   A ....... Matrix                                                      #
#   lrsas ... List of regular semi-algebraic systems                      #
#   R ....... Polynomial ring                                             #
#                                                                         #
# OUTPUT                                                                  #
#     True if all entries in A vanish everywhere over lrsas, false        #
#     otherwise                                                           #
#                                                                         #
# EXAMPLE                                                                 #
#                                                                         #
# ======================================================================= #
# ======================================================================= #
isZeroMatrixOverLRSAS := proc(A::~Matrix, lrsas::TRDlrsas, R::TRDring, $) :: truefalse;

    local n::posint, m::posint, i::posint, j::posint;

    n, m := LA:-Dimension(A);

    for i to n do
        for j to m do
            if not isZeroOverLRSAS(A[i,j], lrsas, R) then
                return false;
            end if;
        end do;
    end do;

    return true;

end proc;