test := module()

    export ModuleApply;

    local test1,
          test2,
          test3,
          test4;

    uses ParametricMatrixTools,
         RegularChains,
         LinearAlgebra;

    ModuleApply := proc($)

        local passCount, failCount, test, testList; 

        testList := ['test1', 'test2', 'test3', 'test4'];

        printf("Testing matrixContainsParameters\n");

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

        local R, A, result, correct;

        R := PolynomialRing([x, a, b]):
        
        A := Matrix(4, 5):

        try
            result := ParametricMatrixTools:-matrixContainsParameters(A, x, R);
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


    test2 := proc($)
    
        local R, A, result, correct;
    
        R := PolynomialRing([x, a, b]):
        
        A := RandomMatrix(4) + x*RandomMatrix(4):
    
        try
            result := ParametricMatrixTools:-matrixContainsParameters(A, x, R);
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
    
        local R, A, result, correct;
    
        R := PolynomialRing([x, a, b]):
        
        A := a*RandomMatrix(4) + b*RandomMatrix(4):
    
        try
            result := ParametricMatrixTools:-matrixContainsParameters(A, x, R);
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
    
        local R, A, result, correct;
    
        R := PolynomialRing([x, a, b]):
        
        A := Matrix(4, 5);
        A[2, 2] := x;
        A[3, 2] := a+4*b-2*x;
    
        try
            result := ParametricMatrixTools:-matrixContainsParameters(A, x, R);
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