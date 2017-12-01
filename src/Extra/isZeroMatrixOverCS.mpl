# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isZeroMatrixOverCS.mpl                                                  #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 9/2017                                                 #
#                                                                         #
# Determine if all entries of a matrix of polynomials vanish at all       #
# points in the zero set of a constructible set.                          #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#     True if all entries in A vanish everywhere over cs.                 #
#                                                                         #
# EXAMPLE                                                                 #
#                                                                         #
# ======================================================================= #
# ======================================================================= #
isZeroMatrixOverCS := proc(A::~Matrix, cs::TRDcs, R::TRDring, $) :: truefalse;

    local n::posint, m::posint, i::posint, j::posint;

    n, m := LA:-Dimension(A);

    for i to n do
        for j to m do
            if not isZeroOverCS(A[i,j], cs, R) then
                return false;
            end if;
        end do;
    end do;

    return true;

end proc;