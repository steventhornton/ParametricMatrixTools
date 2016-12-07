# ======================================================================= #
# ======================================================================= #
#                                                                         #
# TRDdifference_intersect_rs_p.mpl                                        #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 7/2016                                                 #
#                                                                         #
# Compute the difference and intersection of the zero set of a regular    #
# system and the zero set of a polynomial.                                #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   TRDdifference_intersect_rs_p(rs, p, R)                                #
#                                                                         #
# INPUT                                                                   #
#   rs ... Regular system                                                 #
#   p .... Polynomial                                                     #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   Two constructible sets cs_nz, cs_z are returned such that p vanishes  #
#   at all points in the zero set of cs_z and p vanishes nowhere in the   #
#   zero set of cs_nz.                                                    #
#                                                                         #
# EXAMPLE                                                                 #
#   > R := PolynomialRing([a, b]):                                        #
#   > p := a - 1:                                                         #
#   > rs := RegularSystem(Chain([b-1], Empty(R), R), R):                  #
#   > cs_nz, cs_z := TRDdifference_intersect_rs_p(rs, p, R):              #
#   > Display(cs_nz, R);                                                  #
#         { b + 1 = 0                                                     #
#         { a - 1 <> 0                                                    #
#   > Display(cs_z, R);                                                   #
#         { a - 1 = 0                                                     #
#         { b - 1 = 0                                                     #
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
TRDdifference_intersect_rs_p := proc(rs::TRDrs, p::depends(polyInRing(R)), R::TRDring, $) :: TRDcs, TRDcs;
    
    local cs_p::TRDcs, cs_rs::TRDcs;
    
    cs_p := RC_CST:-GeneralConstruct([p], [], R);
    cs_rs := RC_CST:-ConstructibleSet([rs], R);
    
    return RC:-TRDdifference_intersect_cs_cs(cs_rs, cs_p, R);
    
end proc;