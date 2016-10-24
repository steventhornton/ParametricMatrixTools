test := module()

    export ModuleApply;

    local test1,
          test2,
          test3,
          test4,
          test5;

    uses ParametricMatrixTools;

    ModuleApply := proc($)
        printf("Testing isConstantMatrix\n");
        printf("\tTest 1: ...");
        test1();
        printf("\tTest 2: ...");
        test2();
        printf("\tTest 3: ...");
        test3();
        printf("\tTest 4: ...");
        test4();
        printf("\tTest 5: ...");
        test5();
        printf("\n");
    end proc;


    test1 := proc($)

        local A, result;

        A := Matrix(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9]);

        result := ParametricMatrixTools:-isConstantMatrix(A);
        
        printf("\b\b\b");
        if result then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;

    end proc;
    
    
    test2 := proc($)
    
        local A, result;
    
        A := Matrix(3, 3, [1, 2, 3, 4, 5, 6, x, y, z]):
    
        result := ParametricMatrixTools:-isConstantMatrix(A);
        
        printf("\b\b\b");
        if not result then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
    
    end proc;
    
    
    test3 := proc($)
    
        local A, result;
    
        A := Matrix(3, 3, [1, 2, 3, 4, 5, 6, x, y, z]):
    
        result := ParametricMatrixTools:-isConstantMatrix(A, x);
        
        printf("\b\b\b");
        if not result then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
    
    end proc;
    
    
    test4 := proc($)
    
        local A, result;
    
        A := Matrix(3, 4, [x, 2, 3, 4, x, 6, 7, 8, x, 10, 11, 12]):
    
        result := ParametricMatrixTools:-isConstantMatrix(A);
        
        printf("\b\b\b");
        if not result then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
    
    end proc;
    
    
    test5 := proc($)
    
        local A, result;
    
        A := Matrix(3, 4, [x, 2, 3, 4, x, 6, 7, 8, x, 10, 11, 12]):
    
        result := ParametricMatrixTools:-isConstantMatrix(A, x);
        
        printf("\b\b\b");
        if result then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
    
    end proc;


end module: