# ----------------------------------------------------------------------- #
# ListUnion
#
# Computes the union of all elements in a list or set of  constructible sets.
#
# INPUT:
#   lcs ... List or set of constructible sets
#   R ..... Polynomial ring
#
# OUTPUT:
#   A constructible set representing the union of all input sets
# ----------------------------------------------------------------------- #
ListUnion := proc(lcs::{list(TRDcs), set(TRDcs)}, R::TRDring, $) :: TRDcs;

    local n::nonnegint,
          i::posint,
          cs::TRDcs;

    n := nops(lcs);

    if n = 0 then
        return NULL;
    end if;
    if n = 1 then
        return op(lcs);
    end if;

    cs := lcs[1];

    for i from 2 to n do

        cs := RC_CST:-Union(cs, lcs[i], R);

    end do;

    return cs;

end proc: