# ======================================================================= #
# ======================================================================= #
#                                                                         #
# minimal_polynomial_snf.mpl                                              #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 24/2017                                                #
#                                                                         #
# Computes a complete case discussion for the minimal polynomial of a     #
# matrix where the entries are multivariate polynomials whose             #
# indeterminants are regarded as parameters. Computation is done over     #
# polynomial equality and inequation constraints on the parameters.       #
# The minimal polynomial is computed by first computing the Smith form    #
# of A - vI and extracting the minimal polynomial from the SNF. The       #
# algorithm used to compute the SNF is the default for the                #
# ComprehensiveSmithForm routine.                                         #
#                                                                         #
# INPUT                                                                   #
#   A .... Matrix                                                         #
#   v .... Variable                                                       #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [pmin, cs_i]                                                      #
#   where pmin is the minimal polynomial of A for all parameter values    #
#   that satisfy the equations and inequations of cs_i. Together, all the #
#   constructible sets cs_i (for all values of i) represent partition of  #
#   the input constructible set.                                          #
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
comprehensive_minimal_polynomial_snf := module()

    export ModuleApply;

    local
        implementation;
    
    ModuleApply := proc(A::Matrix, v::name, cs::TRDcs, R::TRDring, $)
        return implementation(A, v, cs, R);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# ----------------------------------------------------------------------- #
# implementation                                                          #
#                                                                         #
# Compute the minimal polynomial of a parametric matrix.                  #
#                                                                         #
# INPUT/OUTPUT                                                            #
#    Same as minimal_polynomial_snf                                       #
# ----------------------------------------------------------------------- #
implementation := proc(A::Matrix, v::name, cs::TRDcs, R::TRDring, $)
    
    local n :: posint,
          Ax :: Matrix,
          Rx :: TRDring,
          SList :: list([Matrix, TRDcs]),
          out :: list([polynom, TRDcs]);
    
    out := [];
    
    # Get the size of A
    n := LA:-RowDimension(A);
    
    # Matrix to compute the Smith Form of
    Ax := A - v*LA:-IdentityMatrix(n);
    
    Rx := RC:-PolynomialRing([v, op(R['variables'])]);
    
    # Compute the Smith Normal Form
    SList := ComprehensiveSmithForm(Ax, v, cs, Rx, 'outputType'='CS');
    
    # Extract the minimal polynomial from each matrix in SList
    out := map(x -> [x[1][n,n], x[2]], SList);
    
    return out;
    
end proc;

end module;