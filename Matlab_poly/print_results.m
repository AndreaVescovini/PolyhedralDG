res = fopen('results.txt', 'wt');
fprintf(res, 'EXPONENTIAL - STRUCTURED GRID\n');

uex = @(p) exp(p(1,:).*p(2,:).*p(3,:));
f = @(p) -uex(p).*((p(1,:).*p(2,:)).^2+(p(1,:).*p(3,:)).^2+(p(2,:).*p(3,:)).^2);
uex_grad = @(p) [uex(p).*p(2,:).*p(3,:); uex(p).*p(1,:).*p(3,:); uex(p).*p(1,:).*p(2,:)];

% uex = @(p) p(1,:).*p(2,:);
% f = @(p) 0;
% uex_grad = @(p) [p(2,:); p(1,:); 0];

for N = 1:3
    [err_L2, err_H10, r_L2, r_H10] = hconvergence_ratio(uex, f, uex_grad, N);
    fprintf(res,'\nN = %d\n', N);
    fprintf(res, 'errl2: %.4e %.4e %.4e\n', err_L2);
    fprintf(res, 'errH10: %.4e %.4e %.4e\n', err_H10);
    fprintf(res, 'rl2: %f %f\n', r_L2);
    fprintf(res, 'rH10: %f %f\n', r_H10);
end

fclose(res);