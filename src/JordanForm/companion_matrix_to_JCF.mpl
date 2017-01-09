# ======================================================================= #
# ======================================================================= #
#                                                                         #
# companion_matrix_to_JCF.mpl                                             #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 22/2016                                                #
#                                                                         #
# Compute a full discussion for the Jordan canonical form (JCF) of a      #
# Frobenius companion matrix of a parametric univariate polynomial.       #
#                                                                         #
# INPUT                                                                   #
#   C ... Frobenius companion matrix                                      #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [J_i, cs_i]                                                       #
#   such that J_i is the JCF of the input matrix C for all values of the  #
#   parameters in the zero set of cs_i. The set of constructible sets     #
#   {cs_i, for all i} forms a partition of the input constructible set.   #
#                                                                         #
# EXAMPLE                                                                 #
#   > p := (x + 1) * (x^2 + x + 1) * (x + a):                             #
#   > C := CompanionMatrix(p, x):                                         #
#   > R := PolynomialRing([a]):                                           #
#   > cs := GeneralConstruct([], [], R):                                  #
#   > result := companion_matrix_to_JCF(C, cs, R):                        #
#   > convert(result[1][1], radical), Display(result[1][2], R);           #
#         [-1, 0, 0, 0]                                                   #
#         [0, -1/2+((1/2)*I)*sqrt(3), 0, 0]     { a-1 <> 0                #
#         [0, 0, -1/2-((1/2)*I)*sqrt(3), 0]     { a^2 - a + 1 <> 0        #
#         [0, 0, 0, -a]                                                   #
#   > convert(result[2][1], radical), Display(result[2][2], R);           #
#         [-1/2+((1/2)*I)*sqrt(3), 0, 0, 0]                               #
#         [0, -1/2-((1/2)*I)*sqrt(3), 0, 0]     { a - 1 = 0               #
#         [0, 0, -1, 1]                                                   #
#         [0, 0, 0, -1]                                                   #
#   > convert(result[3][1], radical), Display(result[3][2], R);           #
#         [-1,  0,   0,  0]                                               #
#         [ 0, a-1,  0,  0]     { a^2 - a + 1 = 0                         #
#         [ 0,  0,  -a,  1]                                               #
#         [ 0,  0,   0, -a]                                               #
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
companion_matrix_to_JCF := module()

    export ModuleApply;

    local init,
          implementation,
          constructJCF;

    ModuleApply := proc()
        return init(args);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# ----------------------------------------------------------------------- #
# init                                                                    #
#                                                                         #
# Check the input for errors.                                             #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as comprehensive_FNF_to_JCF                                      #
# ----------------------------------------------------------------------- #
init := proc(A::Matrix(square), cs::TRDcs, R::TRDring, $)
    
    # TO DO
    # Check that A is a matrix of polynomials in R
    
    return implementation(A, cs, R);
    
end proc;


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Compute the full case discusion for the Jordan canonical form of a      #
# companion matrix.                                                       #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as comprehensive_FNF_to_JCF                                      #
# ----------------------------------------------------------------------- #
implementation := proc(A, cs, R, $)
    
    local p :: polynom,
          p_sqr_free_list,
          p_sqr_free,
          q,
          es,
          result,
          R2 :: TRDring,
          x;
    
    # Compute the square-free factorization of the characteristic 
    # polynomial of A
    p := LA:-CharacteristicPolynomial(A, 'x');
    R2 := RC:-PolynomialRing(['x', op(R['variables'])]);
    p_sqr_free_list := ComprehensiveSquareFreeFactorization(p, 'x', cs, R2, 'outputType'='CS');
    
    # Build the JCF of A for each case returned by 
    # ComprehensiveSquareFreeFactorization
    result := [];
    for p_sqr_free in p_sqr_free_list do
        q, es := op(p_sqr_free);
        result := [op(result), [constructJCF(q, 'x'), es]];
    end do;
    
    return result;
    
end proc;


# ----------------------------------------------------------------------- #
# constructJCF                                                            #
#                                                                         #
# Given the square-free decomposition of the characteristic polynomial    #
# of a companion matrix, build it's JCF based on the square-free          #
# decomposition.                                                          #
#                                                                         #
# INPUT                                                                   #
#   sqr_free_list ... A list with elements of the form:                   #
#                         [p_i, m_i]                                      #
#                     where p_i is a polynomial and m_i is a positive     #
#                     integer such that prod((p_i)^(m_i)) is the          #
#                     square-free decomposition in the variable v of some #
#                     parametric univariate polynomial given some         #
#                     contraints on the parameters. m_i <> m_j for i <> j #
#   v ............... A variable                                          #
#                                                                         #
# OUTPUT                                                                  #
#   A Jordan block matrix                                                 #
# ----------------------------------------------------------------------- #
constructJCF := proc(sqr_free_list, v)
    
    local blocks, task, p, m, term, q, n, i;
    
    blocks := [];
    
    for task in sqr_free_list do
        p, m := op(task);
        
        p := sqrfree(factor(p), v)[2];
        
        for term in p do
            q, n := op(term);
            
            ASSERT(n = 1);
            
            if degree(q, v) > 1 then
                for i to degree(q, v) do
                    blocks := [op(blocks), [RootOf(q, v, 'index'=i), m]];
                end do;
            else
                blocks := [op(blocks), [RootOf(q, v), m]];
            end if;
            
        end do;
        
    end do;
    return LA:-JordanBlockMatrix(blocks);
    
end proc;

end module;