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
        printf("Testing ParametricSmithForm_Minors\n");
        
        # determinantDivisor
        printf("\tTesting determinantDivisor\n")
        printf("\t\tTest 1: ...");
        test1();
        printf("\t\tTest 2: ...");
        test2();
        printf("\t\tTest 3: ...");
        test3();
        
        # getMinor
        printf("\tTesting getMinor\n");
        printf("\t\tTest 4: ...");
        test1();
        printf("\t\tTest 5: ...");
        test2();
        printf("\t\tTest 6: ...");
        test3();
        printf("\n");
    end proc;
    
    
    
end module: