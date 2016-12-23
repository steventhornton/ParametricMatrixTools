test := module()

    export ModuleApply;

    local test1, test2, test3, test4, test5,
          test6, test7, test8, test9, test10;

    uses RegularChains, 
         RegularChains:-ConstructibleSetTools,
         RegularChains:-ChainTools,
         ParametricMatrixTools;

    ModuleApply := proc($)

        local passCount, failCount, test, testList; 

        testList := ['test1', 'test2', 'test3', 'test4', 'test5', 
                     'test6', 'test7', 'test8', 'test9', 'test10'];

        printf("Testing TRDis_disjoint_lcs\n");

        passCount, failCount := 0, 0;

        for test in testList do
            printf("\t%a: ...", test);
            if test() then
                passCount := passCount + 1
            else
                failCount := failCount + 1;
            end if;
        end do;

        printf("\n");

        return passCount, failCount;

    end proc;


    test1 := proc($)

        local R, lcs, result, correct;

        R := PolynomialRing([a]);
        lcs := [GeneralConstruct([], [a-2], R),
                GeneralConstruct([a-2], [], R)];

        try
            result := ParametricMatrixTools:-TRDis_disjoint_lcs(lcs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;

        correct := true;

        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;

    end proc;
    
    
    test2 := proc($)
    
        local R, lcs, result, correct;
    
        R := PolynomialRing([a]);
        lcs := [GeneralConstruct([], [a+2, a^3+2], R),
                GeneralConstruct([a+2], [], R),
                GeneralConstruct([a^3+2], [], R)];
    
        try
            result := ParametricMatrixTools:-TRDis_disjoint_lcs(lcs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := true;
    
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test3 := proc($)
    
        local R, lcs, result, correct;
    
        R := PolynomialRing([a,b]);
        lcs := [GeneralConstruct([], [a, a-1, a-b, a*b-b+2, b-1, b+1], R),
                GeneralConstruct([b*a-b+2], [b, b^2-b+2], R),
                GeneralConstruct([a-b], [b-1, b^2-b+2], R),
                GeneralConstruct([a], [b, b-2, b-1], R),
                GeneralConstruct([a-1], [b-1], R),
                GeneralConstruct([b-1], [a, a+1, a-1], R),
                GeneralConstruct([b+1], [a, a-3, a+1, a-1], R),
                GeneralConstruct([b*a-b+2, b^2-b+2], [b], R),
                GeneralConstruct([a-1, b-1], [], R),
                GeneralConstruct([a, b-1], [], R)];
    
        try
            result := ParametricMatrixTools:-TRDis_disjoint_lcs(lcs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := true;
    
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test4 := proc($)
    
        local R, lcs, result, correct;
    
        R := PolynomialRing([a, b]);
        lcs := [GeneralConstruct([], [a^2-b], R),
                GeneralConstruct([a^2-b], [b], R),
                GeneralConstruct([a, b], [b], R)];
    
        try
            result := ParametricMatrixTools:-TRDis_disjoint_lcs(lcs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := true;
    
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test5 := proc($)
    
        local R, lcs, result, correct;
    
        R := PolynomialRing([a, b, c, d]);
        lcs := [GeneralConstruct([a, b, c, d], [], R),
                GeneralConstruct([], [d*a-c*b, d], R),
                GeneralConstruct([d], [b, c], R),
                GeneralConstruct([b, c, d], [a], R),
                GeneralConstruct([a, c, d], [b], R),
                GeneralConstruct([c, d], [a, b], R),
                GeneralConstruct([b, d], [a, c], R),
                GeneralConstruct([a, b, d], [c], R),
                GeneralConstruct([a, c], [d], R),
                GeneralConstruct([d*a-c*b], [b, c, d], R),
                GeneralConstruct([a, b], [c, d], R)];
    
        try
            result := ParametricMatrixTools:-TRDis_disjoint_lcs(lcs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := true;
    
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test6 := proc($)
    
        local R, lcs, result, correct;
    
        R := PolynomialRing([a, b, c, d]);
        lcs := [GeneralConstruct([a, b, c, d], [], R),
                GeneralConstruct([], [d*a-c*b, d], R),
                GeneralConstruct([d], [b, c], R),
                GeneralConstruct([b, c, d], [a], R),
                GeneralConstruct([a, c, d], [b], R),
                GeneralConstruct([c], [a, b], R),
                GeneralConstruct([b, d], [a, c], R),
                GeneralConstruct([a, b, d], [c], R),
                GeneralConstruct([a, c], [d], R),
                GeneralConstruct([d*a-c*b], [b, c, d], R),
                GeneralConstruct([a, b], [c, d], R)];
    
        try
            result := ParametricMatrixTools:-TRDis_disjoint_lcs(lcs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := false;
    
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test7 := proc($)
    
        local R, lcs, result, correct;
    
        R := PolynomialRing([a, b, c]);
        lcs := [GeneralConstruct([a, b, c], [], R),
                GeneralConstruct([a, c], [b], R),
                GeneralConstruct([a+b], [c], R),
                GeneralConstruct([a+b, c], [b], R),
                GeneralConstruct([], [a+b, c], R),
                GeneralConstruct([c], [a, a+b, b], R),
                GeneralConstruct([b, c], [a], R)];
    
        try
            result := ParametricMatrixTools:-TRDis_disjoint_lcs(lcs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := true;
    
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test8 := proc($)
    
        local R, lcs, result, correct;
    
        R := PolynomialRing([a, b]);
        lcs := [GeneralConstruct([a,b], [], R),
                GeneralConstruct([], [a+b,a-b,b], R),
                GeneralConstruct([b], [a], R),
                GeneralConstruct([a+b], [b], R),
                GeneralConstruct([a-b], [b], R)];
    
        try
            result := ParametricMatrixTools:-TRDis_disjoint_lcs(lcs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := true;
    
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test9 := proc($)
    
        local R, lcs, result, correct;
    
        R := PolynomialRing([a, b]);
        lcs := [GeneralConstruct([a+1,b], [], R),
                GeneralConstruct([], [b], R),
                GeneralConstruct([b], [a+1], R)];
    
        try
            result := ParametricMatrixTools:-TRDis_disjoint_lcs(lcs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := true;
    
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;
    
    
    test10 := proc($)
    
        local R, lcs, result, correct;
    
        R := PolynomialRing([a, b, c, d]);
        lcs := [GeneralConstruct([], [a^2*d-a*b*c-a*c*d+b*c^2+b^2-2*b*d+d^2, d], R),
                GeneralConstruct([d], [c*a-c^2-b, b], R),
                GeneralConstruct([a, b, c, d], [], R),
                GeneralConstruct([b, d], [a-c], R),
                GeneralConstruct([c*a-c^2-b, d], [b,c], R),
                GeneralConstruct([a-c, b, d], [c], R),
                GeneralConstruct([d*a^2+(-b*c-c*d)*a+b^2+(c^2-2*d)*b+d^2], [b-d,d], R),
                GeneralConstruct([a-c,b-d], [d], R)];
    
        try
            result := ParametricMatrixTools:-TRDis_disjoint_lcs(lcs, R);
        catch:
            printf("\b\b\bFAIL: Error\n");
            return false;
        end try;
    
        correct := true;
    
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
            return true;
        else
            printf("FAIL: Incorrect result\n");
            return false;
        end if;
    
    end proc;



end module: