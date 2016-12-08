test := module()

    export ModuleApply;

    local test1,
          test2;

    uses ParametricMatrixTools, 
         RegularChains;

    ModuleApply := proc($)
    
        local passCount, failCount, test, testList; 
        
        testList := ['test1','test2'];
        
        printf("Testing isGreatestVariable\n");
        
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

        local R, result, correct;

        R := PolynomialRing([x, a, b]);
        
        try
            result := ParametricMatrixTools:-isGreatestVariable(x, R);
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
    
        local R, result, correct;
    
        R := PolynomialRing([x, a, b]);
        
        try
            result := ParametricMatrixTools:-isGreatestVariable(b, R);
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