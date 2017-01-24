# ======================================================================= #
# ======================================================================= #
#                                                                         #
# FNF_to_JCF.mpl                                                          #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
#                Under the supervision of                                 #
#                Robert M. Corless & Marc Moreno Maza                     #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Jan. 24/2017                                                #
#                                                                         #
# Compute the Jordan form of a matrix in Frobenius form.                  #
#                                                                         #
# INPUT                                                                   #
#   F .... Matrix in Frobenius form                                       #
#   cs ... Constructible set                                              #
#   R .... Polynomial ring                                                #
#                                                                         #
# OUTPUT                                                                  #
#   A list with elements of the form                                      #
#       [J, cs]                                                           #
#   where J is the Jordan form of F for all points in the zero set of cs. #
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
FNF_to_JCF := proc(F::Matrix, cs, R::TRDring, $)

    local charPolyList,
          N,
          R2,
          out,
          tasks,
          JList,
          cs_i,
          i,
          JCFList,
          J2,
          cs_2,
          j;
          
    
    # charPolyList[1] = minPoly(F)
    charPolyList := getCharPolys(F, 'v');
    N := nops(charPolyList);
    
    R2 := RC:-PolynomialRing(['v', op(R['variables'])]);
    
    out := [];
    
    # A list of the form with elements of the form [JList, cs, i]
    # where i is the index of the characteristic polynomial that needs to
    # be computed next
    # JList is a list of matrices
    tasks := poly_to_JCF(numer(charPolyList[1]), 'v', cs, R2);
    
    if N > 1 then
        tasks := map(x -> [[x[1]], x[2], 2], tasks);
    else
        return tasks;
    end if;
    
    while nops(tasks) > 0 do
        
        JList, cs_i, i := op(tasks[1]);
        tasks := tasks[2..-1];
        
        # Compute JCF
        JCFList := poly_to_JCF(numer(charPolyList[i]), 'v', cs_i, R2);
        
        # Add to either output list or tasks
        if i = N then
            for j to nops(JCFList) do
                J2, cs_2 := op(JCFList[j]);
                out := [op(out), [[op(JList), J2], cs_2]];
            end do;
        else
            for j to nops(JCFList) do
                J2, cs_2 := op(JCFList[j]);
                tasks := [op(tasks), [[op(JList), J2], cs_2, i+1]];
            end do;
        end if;
        
    end do;
    
    return map(x -> [LA:-DiagonalMatrix(x[1]), x[2]], out);

end proc;