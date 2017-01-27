# ======================================================================= #
# ======================================================================= #
#                                                                         #
# comprehensive_frobenius_form_snf_minors.mpl                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 9/2017                                                 #
#                                                                         #
# Computes the Frobenius (rational) normal form of a matrix where the     #
# entries are multivariate polynomials by computing the Smith normal      #
# form. Computation is done modulo a regular system or constructible set. #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements in one of the following form:                    #
#       [F, cs]                                                           #
#   Where F is the Frobenius normal form of A for all parameter values    #
#   that satisfy the equations and inequations of cs. Together, all the   #
#   constructible sets form a partition of the input constructible set.   #
#                                                                         #
# ASSUMPTIONS                                                             #
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
comprehensive_frobenius_form_snf_minors := module()

    export ModuleApply;

    local
        implementation,
        SNF_to_FNF;

    ModuleApply := proc(A::~Matrix, cs::TRDcs, R::TRDring, $)
        return implementation(A, cs, R);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Compute the Frobenius Form given a constructible set.                   #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as comprehensive_frobenius_form_snf_minors                       #
# ----------------------------------------------------------------------- #
implementation := proc(A::~Matrix, cs::TRDcs, R::TRDring, $)
    
    local n::posint,
          v::name,
          Ax::Matrix,
          Rx::TRDring,
          SList;
    
    n := LA:-RowDimension(A);
    
    # Matrix to compute the Smith Form of
    Ax := A - v*LA:-IdentityMatrix(n);
    
    Rx := RC:-PolynomialRing([v, op(R['variables'])]);
    
    # Compute the Smith form of Ax
    SList := ComprehensiveSmithForm(Ax, v, cs, Rx, 'outputMatrices'='S', 'outputType'='CS', 'algorithm'='minors');
    
    # Convert matrices in SNF to matrices in Frobenius Form
    return map(x -> [SNF_to_FNF(x[1], v), x[2]], SList);
    
end proc;


# ----------------------------------------------------------------------- #
# SNF_to_FNF                                                              #
#                                                                         #
# Convert a matrix in Smith normal form to a matrix in Frobenius form.    #
#                                                                         #
# INPUT                                                                   #
#   S ... Matrix                                                          #
#   v ... Variable                                                        #
#                                                                         #
# OUTPUT                                                                  #
#    The companion matrices of the polynomials along the diagonal (in     #
#    reverse order) along the diagonal of a matrix (i.e. the Frobenius    #
#    Form).                                                               #
# ----------------------------------------------------------------------- #
SNF_to_FNF := proc(S::Matrix, v::name, $) :: Matrix;

    local CList::list(Matrix) := [],
          n::posint,
          i::posint,
          C::Matrix;

    n := LA:-RowDimension(S);

    for i from n to 1 by -1 do
        if S[i,i] <> 1 and S[i,i] <> 0 then
            C := [LA:-CompanionMatrix(S[i,i], v)][1];
            CList := [op(CList), C];
        end if;
    end do;

    # Put all matrices into a diagonal matrix and return
    return LA:-DiagonalMatrix(CList);

end proc;

end module;