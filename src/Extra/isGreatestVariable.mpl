# ======================================================================= #
# ======================================================================= #
#                                                                         #
# isGreatestVariable.mpl                                                  #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# Determine if a variable is the largest variable in a polynomial ring.   #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   isGreatestVariable(v, R)                                              #
#                                                                         #
# INPUT                                                                   #
#   v ... variable                                                        #
#   R ... Polynomial ring                                                 #
#                                                                         #
# OUTPUT                                                                  #
#   True if v is the largest variable of R, false otherwise.              #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([x, a, b]):                                     #
#   > isGreatestVariable(x, R);                                           #
#         true                                                            #
#   > isGreatestVariable(b, R);                                           #
#         false                                                           #
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
isGreatestVariable := proc(v::name, R::TRDring, $) :: truefalse;
    
    # Ensure v is a variable of R
    ASSERT(evalb(v in R['variables']), "v is not a variable of R");
    
    return evalb(R['variables'][1] = v);
    
end proc;