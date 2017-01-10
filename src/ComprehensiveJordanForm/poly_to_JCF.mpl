# ======================================================================= #
# ======================================================================= #
#                                                                         #
# poly_to_JCF.mpl                                                         #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 10/2017                                                #
#                                                                         #
# Compute a full discussion for the Jordan canonical form (JCF) of a      #
# Frobenius companion matrix of a parametric univariate polynomial.       #
#                                                                         #
# INPUT                                                                   #
#   p .... parametric univariate polynomial                               #
#   v .... Variable                                                       #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [J_i, cs_i]                                                       #
#   such that J_i is the JCF of the Frobenius companion matrix of the     #
#   input polynomial for all parameter values in the zer set of cs_i. The #
#   set of constructible sets {cs_i, for all i} forms a partition of the  #
#   input constructible set.                                              #
#                                                                         #
# EXAMPLE                                                                 #
#   > p := (x + 1) * (x^2 + x + 1) * (x + a):                             #
#   > R := PolynomialRing([x, a]):                                        #
#   > cs := GeneralConstruct([], [], R):                                  #
#   > result := poly_to_JCF(p, x, cs, R):                                 #
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
poly_to_JCF := module()

    export ModuleApply;

    local implementation,
          constructJCF;

    ModuleApply := proc()
        return implementation(args);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Compute the full case discusion for the Jordan canonical form of the    #
# companion matrix associated with a polynomial.                          #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as poly_to_JCF                                                   #
# ----------------------------------------------------------------------- #
implementation := proc(p::depends(polyInRing(R)), v::name, cs, R::TRDring, $)
    
    local p_sqr_free_list,
          p_sqr_free,
          q,
          es,
          result;
    
    # Compute the square-free factorization of p
    p_sqr_free_list := ComprehensiveSquareFreeFactorization(p, v, cs, R, 'outputType'='CS');
    
    # Build the JCF of the companion matrix associated with p for each case 
    # returned by ComprehensiveSquareFreeFactorization
    result := [];
    for p_sqr_free in p_sqr_free_list do
        q, es := op(p_sqr_free);
        result := [op(result), [constructJCF(q, v), es]];
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