test := module()

    export ModuleApply;

    local test1,
          test2,
          test3,
          test4;
    
    uses ParametricMatrixTools;
    
    ModuleApply := proc($)
        printf("Testing allComb\n");
        printf("\tTest 1: ...");
        test1();
        printf("\tTest 2: ...");
        test2();
        printf("\tTest 3: ...");
        test3();
        printf("\tTest 4: ...");
        test4();
        printf("\n");
    end proc;


    test1 := proc($)

        local l, result, correct;
        
        l := [[1, 2, 3], [a, b, c], [x, y]];
        
        result := ParametricMatrixTools:-allComb(l);
        
        correct := [[1, a, x], [2, a, x], [3, a, x], [1, b, x], [2, b, x], [3, b, x], [1, c, x], [2, c, x], [3, c, x], [1, a, y], [2, a, y], [3, a, y], [1, b, y], [2, b, y], [3, b, y], [1, c, y], [2, c, y], [3, c, y]];
        
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
        
    end proc;
    
    
    test2 := proc($)
    
        local l, result, correct;
        
        l := [[1], [2], [3], [4]];
        
        result := ParametricMatrixTools:-allComb(l);
        
        correct := [[1, 2, 3, 4]];
        
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
        
    end proc;
    
    
    test3 := proc($)
    
        local l, result, correct;
        
        l := [[1, 2, 3, 4]];
        
        result := ParametricMatrixTools:-allComb(l);
        
        correct := [[1], [2], [3], [4]];
        
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
        
    end proc;
    
    
    test4 := proc($)
    
        local l, result, correct;
        
        l := [[1, 2, 3, 4, 5], [x]];
        
        result := ParametricMatrixTools:-allComb(l);
        
        correct := [[1, x], [2, x], [3, x], [4, x], [5, x]];
        
        printf("\b\b\b");
        if evalb(result = correct) then
            printf("Pass\n");
        else
            printf("FAIL\n");
        end if;
        
    end proc;
    
end module: