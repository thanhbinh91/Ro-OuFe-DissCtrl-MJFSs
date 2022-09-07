disp('H-inf Performance')
disp(' '); disp('============ alpha = 0.1, beta = 8 ============'); disp(' ');
[P, X, F, L, topt] = SLPMM(0.1, 8);
disp(' '); disp('============ alpha = 0.1, beta = 9 ============'); disp(' ');
[P, X, F, L, topt] = SLPMM(0.1, 9);
disp(' '); disp('============ alpha = 0.1, beta = 10 ============'); disp(' ');
[P, X, F, L, topt] = SLPMM(0.1, 10);