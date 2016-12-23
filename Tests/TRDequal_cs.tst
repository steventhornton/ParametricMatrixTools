test := module()

    export ModuleApply;

    local test1,
          test2,
          test3;

    uses RegularChains, 
         RegularChains:-ConstructibleSetTools,
         RegularChains:-ChainTools,
         ParametricMatrixTools;

    ModuleApply := proc($)

        local passCount, failCount, test, testList; 

        testList := ['test1', 'test2', 'test3'];

        printf("Testing TRDequal_cs\n");

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

        local R, cs1, cs2, es1, es2, es3, result, correct;
        
        R := PolynomialRing([a, b, c]);
        cs1 := GeneralConstruct([a-1, b-1], [c-1], R);
        
        es1 := GeneralConstruct([a-1], [], R);
        es2 := GeneralConstruct([b-1], [], R);
        es3 := GeneralConstruct([], [c-1], R);
        cs2 := ParametricMatrixTools:-ListIntersection([es1, es2, es3], R);
        
        try
            result := ParametricMatrixTools:-TRDequal_cs(cs1, cs2, R);
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
    
        local R, cs1, cs2, es1, es2, result, correct;
        
        R := PolynomialRing([a, b, c]);
        cs1 := GeneralConstruct([a-1, b-1], [c-1], R);
        
        es1 := GeneralConstruct([a-1], [], R);
        es2 := GeneralConstruct([b-1], [], R);
        cs2 := ParametricMatrixTools:-ListIntersection([es1, es2], R);
        
        try
            result := ParametricMatrixTools:-TRDequal_cs(cs1, cs2, R);
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
    
    
    test3 := proc($)
    
        local R, cs1, cs2, es1, es2, es3, result, correct;
        
        R := PolynomialRing([a, b]);
        es1 := GeneralConstruct([a-1], [b-1], R); 
        es2 := GeneralConstruct([a+1], [b+1], R); 
        cs1 := ParametricMatrixTools:-ListUnion([es1, es2], R); 
        cs2 := ParametricMatrixTools:-ListUnion([es2, es1], R);
        
        try
            result := ParametricMatrixTools:-TRDequal_cs(cs1, cs2, R);
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