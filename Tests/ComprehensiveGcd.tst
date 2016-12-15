test := module()

    export ModuleApply;

    local verifyResult_rs,
          verifyResult_cs,
          verifyPartition_cs,
          verifyCofactors,
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
          test11,
          test12,
          test13,
          test14,
          test15,
          test16,
          test17,
          test18,
          test19,
          test21,
          test22,
          test23,
          test24,
          test25,
          test26,
          test27,
          test28,
          test29,
          test30,
          test31,
          test32,
          test33,
          test34,
          test35,
          test36,
          test37,
          test38,
          test39,
          test40,
          test41,
          test42,
          test43,
          test44,
          test45,
          test46,
          test47,
          test48,
          test49,
          test50,
          test51,
          test52,
          test53,
          test54,
          test55,
          test56,
          test57,
          test58,
          test59,
          test60,
          test61,
          test62,
          test63,
          test64,
          test65;

    uses ParametricMatrixTools, RegularChains, RegularChains:-ConstructibleSetTools, RegularChains:-ChainTools;

    ModuleApply := proc($)

        local passCount, failCount, test, testList, t; 

        testList := ['test1',  'test2',  'test3',  'test4',  'test5', 
                     'test6',  'test7',  'test8',  'test9',  'test10', 
                     'test11', 'test12', 'test13', 'test14', 'test15', 
                     'test16', 'test17', 'test18', 'test19', 
                     'test21', 'test22', 'test23', 'test24', 'test25', 
                     'test26', 'test27', 'test28', 'test29', 'test30', 
                     'test31', 'test32', 'test33', 'test34', 'test35', 
                     'test36', 'test37', 'test38', 'test39', 'test40', 
                     'test41', 'test42', 'test43', 'test44', 'test45', 
                     'test46', 'test47', 'test48', 'test49', 'test50', 
                     'test51', 'test52', 'test53', 'test54', 'test55', 
                     'test56', 'test57', 'test58', 'test59', 'test60', 
                     'test61', 'test62', 'test63', 'test64', 'test65'];

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
    verifyPartition_cs := proc(cs, csList, R)
    
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
    
    
    verifyCofactors := proc(p1, p2, result, R)
        
        local item, g, c1, c2, rs;
        
        for item in result do
            g, c1, c2, rs := op(item);
            
            if not ParametricMatrixTools:-isZeroOverRS(numer(c1*g-p1), rs, R) then
                print(SparsePseudoRemainder(c1*g - p1, RepresentingChain(rs, R), R));
                return false;
            end if;
            if not ParametricMatrixTools:-isZeroOverRS(numer(c2*g-p2), rs, R) then
                print(SparsePseudoRemainder(c2*g - p2, RepresentingChain(rs, R), R));
                return false;
            end if;
            
        end do;
        
        return true;
        
    end proc;
    
    
    verifyResult_rs := proc(p1, p2, g, cs, cs_zero, R)
        
        local csList;
        
        # Verify that p1 and p2 both vanish for all values in the zero set 
        # of cs_zero.
        if not IsEmpty(cs_zero, R) then
            if not ParametricMatrixTools:-isZeroOverCS(p1, cs_zero, R) then
                printf("p1 is not zero over cs_zero\n");
                return false;
            end if;
            if not ParametricMatrixTools:-isZeroOverCS(p2, cs_zero, R) then
                printf("p2 is not zero over cs_zero\n");
                return false;
            end if;
        end if;
        
        # Verify that the regular systems in g, plus cs_zero form a 
        # partition of cs
        csList := ListTools:-Flatten([map(x -> x[-1], g)]);
        csList := map(x -> ConstructibleSet([x], R), csList);
        csList := [op(csList), cs_zero];
        if not verifyPartition_cs(cs, csList, R) then
            printf("Not a partition!");
            return false;
        end if;
        
        # If the cofactors are computed, verify that they are correct
        if nops(g[1]) = 4 then
            if not verifyCofactors(p1, p2, g, R) then
                printf("Cofactors are incorrect");
                return false;
            end if;
        end if;
        
        return true;
        
    end proc;
    
    
    verifyResult_cs := proc(p1, p2, g, cs, cs_zero, R)
    
        local pair, G, lrs, cs_p, rs, rc, r1, r2, correct, csList;
        
        # Verify that p1 and p2 both vanish for all values in the zero set 
        # of cs_zero.
        if not IsEmpty(cs_zero, R) then
            if not ParametricMatrixTools:-isZeroOverCS(p1, cs_zero, R) then
                printf("p1 is not zero over cs_zero\n");
                return false;
            end if;
            if not ParametricMatrixTools:-isZeroOverCS(p2, cs_zero, R) then
                printf("p2 is not zero over cs_zero\n");
                return false;
            end if;
        end if;
        
        # Verify that the constructible sets in g, plus cs_zero form a 
        # partition of cs
        csList := ListTools:-Flatten([map(x -> x[-1], g), cs_zero]);
        if not verifyPartition_cs(cs, csList, R) then
            printf("Not a partition!");
            return false;
        end if;
        
        
        for pair in g do
            G, cs_p := op(pair);
            lrs := RepresentingRegularSystems(cs_p, R);
            for rs in lrs do
                rc := RepresentingChain(rs, R);
                
                # Verify that mG = pq + r for p = p1 and p2 and r is always
                # zero
                r1 := sprem(p1, G, x);
                r2 := sprem(p2, G, x);
                
                if not ParametricMatrixTools:-isZeroOverRS(r1, rs, R) then
                    return false;
                end if;
                if not ParametricMatrixTools:-isZeroOverRS(r2, rs, R) then
                    return false;
                end if;
                
            end do;
        end do;
        
        return true;
    
    end proc;


    # ------------------------------------------------------------------- #
    # TEST 1-5                                                            #
    #                                                                     #
    # p1 = (x+2)*(x+a)                                                    #
    # p2 = (x+2)^2                                                        #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
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
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
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
    
        p1 := (x+2)*(x+a);
        p2 := (x+2)^2;
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
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
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := (x+2)*(x+a);
        p2 := (x+2)^2;
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
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
    
        local p1, p2, R, rs, cs, g, correct, cs_zero;
    
        p1 := (x+2)*(x+a);
        p2 := (x+2)^2;
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 6-10                                                           #
    #                                                                     #
    # p1 = (x-1)*(x+a)                                                    #
    # p2 = a*x^2 + 2                                                      #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test6 := proc($)
    
        local p1, p2, R, rs, cs, g, correct, cs_zero;
        
        p1 := (x-1)*(x+a);
        p2 := a*x^2 + 2;
        
        R := PolynomialRing([x,a]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([],[],R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
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
    
        local p1, p2, R, cs, g, correct, cs_zero;
        
        p1 := (x-1)*(x+a);
        p2 := a*x^2 + 2;
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
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
        
        p1 := (x-1)*(x+a);
        p2 := a*x^2 + 2;
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
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
    
        local p1, p2, R, cs, g, correct, cs_zero;
        
        p1 := (x-1)*(x+a);
        p2 := a*x^2 + 2;
        
        R := PolynomialRing([x,a]);
        cs := GeneralConstruct([],[],R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 11-15                                                          #
    #                                                                     #
    # p1 = (x+a)*(x-b+2)*(x+1)                                            #
    # p2 = (x+b)*(x+1)^2*(x-a*b);                                         #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test11 := proc($)
    
        local p1, p2, R, rs, g, correct, cs_zero, cs;
    
        p1 := (x+a)*(x-b+2)*(x+1);
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test12 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := (x+a)*(x-b+2)*(x+1):
        p2 := (x+b)*(x+1)^2*(x-a*b);
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test13 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := (x+a)*(x-b+2)*(x+1):
        p2 := (x+b)*(x+1)^2*(x-a*b);
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test14 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := (x+a)*(x-b+2)*(x+1):
        p2 := (x+b)*(x+1)^2*(x-a*b);
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test15 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := (x+a)*(x-b+2)*(x+1);
        p2 := (x+b)*(x+1)^2*(x-a*b);
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 16-19                                                          #
    #                                                                     #
    # p1 = (x+2)*(x+a)                                                    #
    # p2 = (x+2)^2                                                        #
    # cs = {a<>3, a<>2}                                                   #
    # ------------------------------------------------------------------- #
    test16 := proc($)
    
        local p1, p2, R, rs, cs, g, correct, cs_zero;
    
        p1 := (x+2)*(x+a);
        p2 := (x+2)^2;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[a-3, a-2], R);
        rs := RepresentingRegularSystems(cs, R)[1];
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test17 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := (x+2)*(x+a);
        p2 := (x+2)^2;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[a-3, a-2], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test18 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := (x+2)*(x+a);
        p2 := (x+2)^2;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[a-3, a-2], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [a-3, a-2], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test19 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := (x+2)*(x+a);
        p2 := (x+2)^2;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[a-3, a-2], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 21-25                                                          #
    #                                                                     #
    # p1 = x + 1                                                          #
    # p2 = a*x + b                                                        #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test21 := proc($)
    
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test22 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := x+a;
        p2 := a*x+b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[],R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test23 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := x+a;
        p2 := a*x+b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[],R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test24 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := x+a;
        p2 := a*x+b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[],R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test25 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := x+a;
        p2 := a*x+b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([],[],R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 26-30                                                          #
    #                                                                     #
    # p1 = a*x + b                                                        #
    # p2 = c*x + d                                                        #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test26 := proc($)
    
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test27 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := a*x+b;
        p2 := c*x+d;
        
        R := PolynomialRing([x,a,b,c,d]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test28 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := a*x+b;
        p2 := c*x+d;
        
        R := PolynomialRing([x,a,b,c,d]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test29 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := a*x+b;
        p2 := c*x+d;
        
        R := PolynomialRing([x,a,b,c,d]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test30 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := a*x+b;
        p2 := c*x+d;
        
        R := PolynomialRing([x,a,b,c,d]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    
    # ------------------------------------------------------------------- #
    # TEST 31-35                                                          #
    #                                                                     #
    # p1 = a + b                                                          #
    # p2 = a*x^2 + c                                                      #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test31 := proc($)
    
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test32 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := a + b;
        p2 := a*x^2 + c;
        
        R := PolynomialRing([x,a,b,c]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test33 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := a + b;
        p2 := a*x^2 + c;
        
        R := PolynomialRing([x,a,b,c]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test34 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := a + b;
        p2 := a*x^2 + c;
        
        R := PolynomialRing([x,a,b,c]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test35 := proc($)
    
        local p1, p2, R, g, correct, cs_zero, cs;
    
        p1 := a + b;
        p2 := a*x^2 + c;
        
        R := PolynomialRing([x,a,b,c]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 36-40                                                          #
    #                                                                     #
    # p1 = 2                                                              #
    # p2 = x + c                                                          #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test36 := proc($)
    
        local p1, p2, R, rs, cs, g, correct, cs_zero;
    
        p1 := 2;
        p2 := x + c;
        
        R := PolynomialRing([x,a,b,c]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test37 := proc($)
    
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test38 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := 2;
        p2 := x + c;
        
        R := PolynomialRing([x,a,b,c]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test39 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := 2;
        p2 := x + c;
        
        R := PolynomialRing([x,a,b,c]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test40 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := 2;
        p2 := x + c;
        
        R := PolynomialRing([x,a,b,c]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 41-45                                                          #
    #                                                                     #
    # p1 = b*x + a                                                        #
    # p2 = a*x + b                                                        #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test41 := proc($)
    
        local p1, p2, R, rs, cs, g, correct, cs_zero;
    
        p1 := b*x + a;
        p2 := a*x + b;
        
        R := PolynomialRing([x,a,b]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test42 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := b*x + a;
        p2 := a*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test43 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := b*x + a;
        p2 := a*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test44 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := b*x + a;
        p2 := a*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test45 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero;
    
        p1 := b*x + a;
        p2 := a*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 46-50                                                          #
    #                                                                     #
    # p1 = 5                                                              #
    # p2 = (a+1)*x + b                                                    #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test46 := proc($)
    
        local p1, p2, R, rs, cs, g, correct, cs_zero, csCorrect;
    
        p1 := 5;
        p2 := (a+1)*x + b;
        
        R := PolynomialRing([x,a,b]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test47 := proc($)
    
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test48 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect;
    
        p1 := 5;
        p2 := (a+1)*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test49 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect;
    
        p1 := 5;
        p2 := (a+1)*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test50 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect;
    
        p1 := 5;
        p2 := (a+1)*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 51-55                                                          #
    #                                                                     #
    # p1 = 0                                                              #
    # p2 = (a+1)*x + b                                                    #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test51 := proc($)
    
        local p1, p2, R, rs, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := 0;
        p2 := (a+1)*x + b;
        
        R := PolynomialRing([x,a,b]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1, b], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test52 := proc($)
    
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
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test53 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := 0;
        p2 := (a+1)*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1, b], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test54 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := 0;
        p2 := (a+1)*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1, b], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test55 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := 0;
        p2 := (a+1)*x + b;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 56-60                                                          #
    #                                                                     #
    # p1 = (a+1)*x + (a+1)                                                #
    # p2 = a+1                                                            #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test56 := proc($)
    
        local p1, p2, R, rs, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := (a+1)*x + (a+1);
        p2 := (a+1);
        
        R := PolynomialRing([x,a,b]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    test57 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := (a+1)*x + (a+1);
        p2 := (a+1);
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test58 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := (a+1)*x + (a+1);
        p2 := (a+1);
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test59 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := (a+1)*x + (a+1);
        p2 := (a+1);
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test60 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := (a+1)*x + (a+1);
        p2 := (a+1);
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    # ------------------------------------------------------------------- #
    # TEST 61-65                                                          #
    #                                                                     #
    # p1 = 2*a + 2                                                        #
    # p2 = a + 1                                                          #
    # cs = {}                                                             #
    # ------------------------------------------------------------------- #
    test61 := proc($)
    
        local p1, p2, R, rs, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := 2*a + 2;
        p2 := a + 1;
        
        R := PolynomialRing([x,a,b]);
        rs := RegularSystem(Empty(R), R);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, rs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    test62 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := 2*a + 2;
        p2 := a + 1;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test63 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := 2*a + 2;
        p2 := a + 1;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test64 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := 2*a + 2;
        p2 := a + 1;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, [], [], R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        cs_zero_correct := GeneralConstruct([a+1], [], R);
        if not equalCS(cs_zero_correct, cs_zero, R) then
            printf("\b\b\bFAIL: Incorrect result\n");
            return false;
        end if;
        
        correct := verifyResult_cs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test65 := proc($)
    
        local p1, p2, R, cs, g, correct, cs_zero, csCorrect, cs_zero_correct;
    
        p1 := 2*a + 2;
        p2 := a + 1;
        
        R := PolynomialRing([x,a,b]);
        cs := GeneralConstruct([], [], R);
        
        try
            g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS', 'cofactors'=true);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
        
        correct := verifyResult_rs(p1, p2, g, cs, cs_zero, R);
        
        printf("\b\b\b");
        if correct then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc


end module: