% for beta = 1.9:-0.02:0
%     alpha = 0;
%     fprintf('\n Solving performance with alpha = %4.3f, beta = %4.3f\n', alpha, beta);
%     [X, X_, P, F, L, copt] = SLPMM(alpha, beta);
%     if copt > 8.01
%         break
%     end
%     save('P','P');
%     save('F','F');
%     save('L','L');
% end
% 
% 
% for beta = 3:-0.02:0
%     alpha = 0.1;
%     fprintf('\n Solving performance with alpha = %4.3f, beta = %4.3f\n', alpha, beta);
%     [X, X_, P, F, L, copt] = SLPMM(alpha, beta);
%     if copt > 8.01
%         break
%     end
%     save('P','P');
%     save('F','F');
%     save('L','L');
% end


for beta = 5:-0.02:0
    alpha = 0.2;
    fprintf('\n Solving performance with alpha = %4.3f, beta = %4.3f\n', alpha, beta);
    [X, X_, P, F, L, copt] = SLPMM(alpha, beta);
    if copt > 8.01
        break
    end
    save('P','P');
    save('F','F');
    save('L','L');
end