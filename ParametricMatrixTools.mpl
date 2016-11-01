# ======================================================================= #
# ======================================================================= #
#                                                                         #
# ParametricMatrixTools.mpl                                               #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Nov. 1/2016                                                 #
#                                                                         #
# A module for computations on parametric matrices.                       #
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
ParametricMatrixTools := module()

    option package;
    
    export isUnder;
    
    local ModuleLoad,
          loadTypes,
          allComb,
          containsParameters,
          InitialByVar,
          isConstantMatrix,
          isEqualOverCS,
          isGreatestVariable,
          isNonZeroOverCS,
          isNonZeroOverRS,
          isZeroMatrix,
          isZeroOverCS,
          isZeroOverRS,
          ListIntersection,
          matrixContainsParameters,
          splitConstructibleSet,
          TailByVar,
          VanishingConstructibleSet;
    
    # This function is run when the package is loaded.
    ModuleLoad := proc()
        kernelopts('opaquemodules' = false);
        loadTypes();
    end proc:


# External Files
$include "types.mpl"

$include "Extra/allComb.mpl"
$include "Extra/containsParameters.mpl"
$include "Extra/InitialByVar.mpl"
$include "Extra/isConstantMatrix.mpl"
$include "Extra/isEqualOverCS.mpl"
$include "Extra/isGreatestVariable.mpl"
$include "Extra/isNonZeroOverCS.mpl"
$include "Extra/isNonZeroOverRS.mpl"
$include "Extra/isUnder.mpl"
$include "Extra/isZeroMatrix.mpl"
$include "Extra/isZeroOverCS.mpl"
$include "Extra/isZeroOverRS.mpl"
$include "Extra/ListIntersection.mpl"
$include "Extra/matrixContainsParameters.mpl"
$include "Extra/splitConstructibleSet.mpl"
$include "Extra/TailByVar.mpl"
$include "Extra/VanishingConstructibleSet.mpl"

end module: