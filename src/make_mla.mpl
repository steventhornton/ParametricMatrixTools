with(LibraryTools):

loc := cat(currentdir(),"/ParametricMatrixTools.mla"):

Create(loc):

read("src/ParametricMatrixTools.mpl"):

Save(ParametricMatrixTools, loc):