test := module()

    export ModuleApply;

    local verify,
          verifyPartition,
          equalCS,
          test1,
          test2,
          test3,
          test4,
          test5,
          test6,
          test7,
          test8,
          test9,
          test10,
          test11;

    uses ParametricMatrixTools, RegularChains, RegularChains:-ConstructibleSetTools, RegularChains:-ChainTools;

    ModuleApply := proc($)

        local passCount, failCount, test, testList, t; 

        testList := ['test1', 'test2', 'test3', 'test4', 'test5', 'test6', 'test7', 'test8', 'test9', 'test10', 'test11'];

        printf("Testing ParametricGcd\n");

        passCount, failCount := 0, 0;

        for test in testList do
            printf("\t%a: ...", test);
            try
                t := test();
            catch:
                printf("\b\b\b");
                printf("Error in test!\n");
                t := false;
            end try;
            if t then
                passCount := passCount + 1
            else
                failCount := failCount + 1;
            end if;
        end do;

        printf("\n");

        return passCount, failCount;

    end proc;
    
    equalCS := proc(cs1, cs2, R)
        
        local csD;
        
        csD := Difference(cs1, cs2, R);
        if not IsEmpty(csD, R) then
            return false;
        end if;
        
        csD := Difference(cs2, cs1, R);
        if not IsEmpty(csD, R) then
            return false;
        end if;
        
        return true;
        
    end proc;
    
    
    # Verify that a list of constructible sets forms a partition of
    # another constructible set.
    verifyPartition := proc(cs, csList, R)
    
        local csU, csListNonEmpty, i, j, csI;
        
        csListNonEmpty := remove(c -> IsEmpty(c, R), csList);
        
        csU := ParametricMatrixTools:-ListUnion(csListNonEmpty, R);
        
        if not equalCS(cs, csU, R) then
            return false;
        end if;
        
        # Verify that each pairwise element of csListNonEmpty is disjoint
        if nops(csListNonEmpty) = 1 then
            return true;
        end if;
        for i to nops(csListNonEmpty)-1 do
            for j from i+1 to nops(csListNonEmpty) do
                csI := Intersection(csListNonEmpty[i], csListNonEmpty[j], R);
                if not IsEmpty(csI, R) then
                    return false;
                end if;
            end do;
        end do;
        
        return true;
        
    end proc;
    
    
    verify := proc(p1, p2, g, cs, cs_zero, R)
    
        local pair, G, lrs, cs_p, rs, rc, p1_s, p2_s, r1, r2, correct, csList;
        
        # Verify that p1 and p2 are zero in cs_zero
        if not IsEmpty(cs_zero, R) then
            # printf('Verify p1 and p2 are zero over cs_zero\n');
            if not ParametricMatrixTools:-isZeroOverCS(p1, cs_zero, R) then
                printf("p1 is not zero over cs_zero\n");
                return false;
            end if;
            if not ParametricMatrixTools:-isZeroOverCS(p2, cs_zero, R) then
                printf("p2 is not zero over cs_zero\n");
                return false;
            end if;
        end if;
        
        csList := ListTools:-Flatten([map(x -> x[2], g), cs_zero]);
        
        if not verifyPartition(cs, csList, R) then
            printf("Not a partition!");
            return false;
        end if;
        
        # Verify that the constructible sets in g, plus cs_zero form a 
        # partition of cs
        correct := true;
        for pair in g do
            G, cs_p := op(pair);
            lrs := RepresentingRegularSystems(cs_p, R);
            for rs in lrs do
                rc := RepresentingChain(rs, R);
                
                p1_s := SparsePseudoRemainder(p1, rc, R);
                p2_s := SparsePseudoRemainder(p2, rc, R);
                
                if ParametricMatrixTools:-isZeroOverRS(G, rs, R) then
                    if not ParametricMatrixTools:-isZeroOverRS(p1_s, rs, R) or
                       ParametricMatrixTools:-isZeroOverRS(p2_s, rs, R) then
                       correct := false;
                    end if;
                else
                    r1 := sprem(p1_s, G, x);
                    r2 := sprem(p2_s, G, x);
                    
                    if not ParametricMatrixTools:-isZeroOverRS(r1, rs, R) then
                        correct := false;
                    end if;
                    if not ParametricMatrixTools:-isZeroOverRS(r2, rs, R) then
                        correct := false;
                    end if;
                    
                end if;
                
            end do;
        end do;
    
        return correct;
    
    end proc;


    test1 := proc($)

        local p1, p2, R, rs, cs, g, correct, cs_zero;

        p1 := (x+2)*(x+a);
        p2 := (x+2)^2;
        
        R := PolynomialRing([x,a]);
        rs := RegularSystem(Empty(R), [], R);
        cs := ConstructibleSet([rs], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;
    
    
    test2 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := (x+2)*(x+a);
        p2 := (x+2)^2;
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test3 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
        
        p1 := (x-1)*(x+a);
        p2 := a*x^2 + 2;
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test4 := proc($)
    
        local p1, p2, R, rs, g, correct, cs_zero, cs;
    
        p1 := (x+a)*(x-b+2)*(x+1):
        p2 := (x+b)*(x+1)^2*(x-a*b);
        
        R := PolynomialRing([x,a,b]);
        rs := RegularSystem(Empty(R), R);
        cs := ConstructibleSet([rs], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test5 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := (x+1)*(x+a);
        p2 := (x+2)^2;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[a-3, a-2], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test6 := proc($)
    
        local p1, p2, R, rs, g, correct, cs_zero, cs;
    
        p1 := x+a;
        p2 := a*x+b;
        
        R := PolynomialRing([x,a,b]);
        rs := RegularSystem(Empty(R), R);
        cs := ConstructibleSet([rs], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test7 := proc($)
    
        local p1, p2, R, rs, g, correct, cs_zero, cs;
    
        p1 := a*x+b;
        p2 := c*x+d;
        
        R := PolynomialRing([x,a,b,c,d]);
        rs := RegularSystem(Empty(R), R);
        cs := ConstructibleSet([rs], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test8 := proc($)
    
        local p1, p2, R, rs, g, correct, cs_zero, cs;
    
        p1 := a + b;
        p2 := a*x^2 + c;
        
        R := PolynomialRing([x,a,b,c]);
        rs := RegularSystem(Empty(R), R);
        cs := ConstructibleSet([rs], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test9 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := 2;
        p2 := x + c;
        
        R := PolynomialRing([x,a,b,c]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test10 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect;
    
        p1 := 5;
        p2 := (a+1)*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        if not IsEmpty(cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        if g[1][1] <> 1 then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        csCorrect := GeneralConstruct([],[],R);
        
        if not equalCS(csCorrect, g[1][2], R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test11 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := 0;
        p2 := (a+1)*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1, b], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verify(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;


end module: