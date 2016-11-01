test := module()

    export ModuleApply;

    local test1,
          test2,
          test3;

    uses ParametricMatrixTools;

    ModuleApply := proc($)
    
        local passCount, failCount, test, testList; 
        
        testList := ['test1','test2', 'test3'];
        
        printf("Testing isConstantMatrix\n");
        
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

        local A, result, correct;

        A := Matrix(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9]);
        
        try
            result := ParametricMatrixTools:-isConstantMatrix(A);
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
    
        local A, result, correct;
    
        A := Matrix(3, 3, [1, 2, 3, 4, 5, 6, x, y, z]):
    
        try
            result := ParametricMatrixTools:-isConstantMatrix(A);
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
    
        local A, result, correct;
    
        A := Matrix(3, 4, [x, 2, 3, 4, x, 6, 7, 8, x, 10, 11, 12]):
    
        try
            result := ParametricMatrixTools:-isConstantMatrix(A);
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