# ======================================================================= #
# ======================================================================= #
#                                                                         #
# TRDequal_cs.mpl                                                         #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Dec. 23/2016                                                #
#                                                                         #
# Determine if two constructible sets are equal.                          #
#                                                                         #
# INPUT                                                                   #
#   cs1 ... Constructible set                                             #
#   cs2 ... Constructible set                                             #
#   R ..... Polynomial ring                                               #
#                                                                         #
# OUTPUT                                                                  #
#   True if cs1 equals cs2, false otherwise.                              #
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
TRDequal_cs := proc(cs1::TRDcs, cs2::TRDcs, R::TRDring, $) :: truefalse;
    
    return RC:-TRDconstructible_set_is_contained(cs1, cs2, R) and RC:-TRDconstructible_set_is_contained(cs2, cs1, R);
    
end proc;