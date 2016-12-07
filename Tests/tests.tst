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
if passCount <> 0 then mP := ceil(log10(passCount+1)) else mP := 1 end if:
if failCount <> 0 then mF := ceil(log10(failCount+1)) else mF := 1 end if:
m := max(mP, mF):
printf(cat("%", m, "d Passed\n"), passCount);
printf(cat("%", m, "d Failed\n"), failCount);