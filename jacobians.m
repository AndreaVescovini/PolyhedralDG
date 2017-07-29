function [Jdet, Jcof] = jacobians(x, y, z, r, s, t)
%jacobians functions that computes for every element the absolute value of
%the determinant of the jacobian and the matrix jcof wrt the simplex (-1,-1,-1),
%(-1,-1,1), (-1,1,-1), (1,-1,-1) (pag 176 of
%Modellistica numerica per problemi differenziali, Quarteroni).

K = size(x,2);
Jdet = zeros(K,1);
Jcof = zeros(3,3,K);

for k = 1:K
%     jac_trasp = zeros(3,3);
%     jac_trasp(1,:) = [r(2:end)-r(1) s(2:end)-s(1) t(2:end)-t(1)]\(x(2:end,k)-x(1,k));
%     jac_trasp(2,:) = [r(2:end)-r(1) s(2:end)-s(1) t(2:end)-t(1)]\(y(2:end,k)-y(1,k));
%     jac_trasp(3,:) = [r(2:end)-r(1) s(2:end)-s(1) t(2:end)-t(1)]\(z(2:end,k)-z(1,k)); % i +1 sono per la traslazione dell'elemento di riferimento
    
    jac_trasp = zeros(3,4);
    jac_trasp(1,:) = [r s t ones(4,1)]\x(:,k);
    jac_trasp(2,:) = [r s t ones(4,1)]\y(:,k);
    jac_trasp(3,:) = [r s t ones(4,1)]\z(:,k); % i +1 sono per la traslazione dell'elemento di riferimento
    
    jac_trasp = jac_trasp(:,1:3);
    
    Jdet(k) = det(jac_trasp);
    Jcof(:,:,k) = Jdet(k)*inv(jac_trasp);
    
end

end

