kernelopts(opaquemodules = false):

printf("# ------------------------------ #\n");
printf("# Running All Tests              #\n");
printf("# ------------------------------ #\n\n");

passCount := 0:
failCount := 0:

# Get all tst files in the Tests directory
testFiles := FileTools:-ListDirectory("Tests/", 'returnonly' = "*.tst"):

# Run each test file
for testFile in testFiles do
    if testFile <> "tests.tst" then
        read(cat("Tests/", testFile)):
        p, f := test():
        passCount := passCount + p:
        failCount := failCount + f:
    end if:
end do:

# Print a summary of the tests
printf("%d Passed\n", passCount);
printf("%d Failed\n", failCount);