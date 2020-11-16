A=[3//5 6//5 3//1 1 0 0; 3//7 4//5 3//2 0 1 0; 7//3 8//5 1 0 0 1];
b=[6//1; 5//1; 7//1];
cj = [3//1 5//1 6//1 0 0 0];
x = ["x1", "x2", "x3", "s1", "s2", "s3"];
xB = ["s1", "s2", "s3"];

include("simplex.jl");
A, b, cj, x, xB, cB, z, zc, ratio, pivColIdx, pivRowIdx = simplexIterator(A, b, cj, x, xB);
