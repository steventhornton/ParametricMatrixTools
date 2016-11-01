test := module()

    export ModuleApply;

    local test1,
          test2,
          test3;

    uses ParametricMatrixTools, 
         RegularChains, 
         RegularChains:-ConstructibleSetTools;

    ModuleApply := proc($)
    
        local passCount, failCount, test, testList; 
        
        testList := ['test1','test2', 'test3'];
        
        printf("Testing isEqualOverCS\n");
        
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

        local R, p1, p2, cs, result, correct;
        
        R := PolynomialRing([x, a]);
        
        p1 := 2*x^2 - 12;
        p2 := (a-1)*x^2 - 4*a;
        
        cs := GeneralConstruct([a-3], [], R):
        
        try
            result := ParametricMatrixTools:-isEqualOverCS(p1, p2, cs, R);
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
    
        local R, p1, p2, cs, result, correct;
        
        R := PolynomialRing([x, a]);
        
        p1 := 2*x^2 - 12;
        p2 := (a-1)*x^2 - 4*a;
        
        cs := GeneralConstruct([], [a-3], R):
        
        try
            result := ParametricMatrixTools:-isEqualOverCS(p1, p2, cs, R);
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
    
        local R, p1, p2, cs, result, correct;
        
        R := PolynomialRing([x, a]);
        
        p1 := 2*x^2 - 12;
        p2 := (a-1)*x^2 - 4*a;
        
        cs := GeneralConstruct([], [], R):
        
        try
            result := ParametricMatrixTools:-isEqualOverCS(p1, p2, cs, R);
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

end module: