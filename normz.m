%Matrix Norm for Matrix A
%L1 Norm
A = magic(10);
B= [1,3,4,5,6,7,2,33,4,55,66,66,77,77,8,32,2,2,2];
l1 = max(sum(abs(A)));
matL1 = norm(A,1);
%L2 Norm
l2 = max(svd(A));
matL2 = norm(A,2);
%Comparison
fprintf('L1 Norm | %g\nMATLAB L1 Norm | %g\n',l1,matL1);
fprintf('L2 Norm | %g\nMATLAB L2 Norm | %g\n',l2,matL2);
%Vector Norm for Vector B
%L1 Norm
l1 = (sum(abs(B)));
matL1 = norm(B,1);
%L2 Norm
l2 = sqrt(sum(abs(B).^2));
matL2 = norm(B,2);
%Comparison
fprintf('\nL1 Norm | %g\nMATLAB L1 Norm | %g\n',l1,matL1);
fprintf('L2 Norm | %g\nMATLAB L2 Norm | %g\n',l2,matL2);