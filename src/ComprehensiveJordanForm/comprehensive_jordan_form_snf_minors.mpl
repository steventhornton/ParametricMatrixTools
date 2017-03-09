# ======================================================================= #
# ======================================================================= #
#                                                                         #
# comprehensive_jordan_form_snf_minors.mpl                                #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 29/2017                                                #
#                                                                         #
# Computes the Jordan canonical form of a matrix where the entries are    #
# mutlivariate polynomials by computing the Smith form of xI - A.         #
# Computation is done modulo a regular system or constructible set.       #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements in one of the following form:                    #
#       [J, cs]                                                           #
#   Where J is the Jordan canonical form of A for all parameter values    #
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
comprehensive_jordan_form_snf_minors := module()

    export ModuleApply;

    local
        implementation,
        FNF_to_JCF,
        poly_to_JCF,
        getCharPolys;

    ModuleApply := proc(A::~Matrix, cs::TRDcs, R::TRDring, $)
        return implementation(A, cs, R);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Compute the Jordan form given a constructible set.                      #
#                                                                         #
# INPUT/OUTPUT                                                            #
#   Same as comprehensive_jordan_form_snf_minors                          #
# ----------------------------------------------------------------------- #
implementation := proc(A::~Matrix, cs::TRDcs, R::TRDring, $)

    local FList;
    
    # Compute the Frobenius form of A
    FList := ComprehensiveFrobeniusForm(A, cs, R, 'outputMatrices'='F', 'outputType'='CS', 'algorithm'='snf_minors');
    
    # For element of FList, compute the Jordan form
    return map(x -> op(FNF_to_JCF(x[1], x[2], R)), FList);

end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# EXTERNAL FILES
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$include <src/ComprehensiveJordanForm/FNF_to_JCF.mpl>
$include <src/ComprehensiveJordanForm/poly_to_JCF.mpl>
$include <src/ComprehensiveJordanForm/getCharPolys.mpl>

end module;