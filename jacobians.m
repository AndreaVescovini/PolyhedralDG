function [J, Jcof, Jdet, trasl] = jacobians(x, y, z, r, s, t)
%jacobians functions that computes for every element the determinant of
%the jacobian and the matrix jcof wrt the simplex (0,0,0),(1,0,0),(0,1,0),(0,0,1)
%(pag 176 of Modellistica numerica per problemi differenziali, Quarteroni).

K = size(x,2);
J = zeros(3,3,K);
Jcof = zeros(3,3,K);
trasl = zeros(3,K);
Jdet = zeros(K,1);

for k = 1:K
%     jac = zeros(3,3);
%     jac(1,:) = [r(2:end)-r(1) s(2:end)-s(1) t(2:end)-t(1)]\(x(2:end,k)-x(1,k));
%     jac(2,:) = [r(2:end)-r(1) s(2:end)-s(1) t(2:end)-t(1)]\(y(2:end,k)-y(1,k));
%     jac(3,:) = [r(2:end)-r(1) s(2:end)-s(1) t(2:end)-t(1)]\(z(2:end,k)-z(1,k)); % i +1 sono per la traslazione dell'elemento di riferimento
    
    Jtemp = zeros(3,4);
    Jtemp(1,:) = [r s t ones(4,1)]\x(:,k);
    Jtemp(2,:) = [r s t ones(4,1)]\y(:,k);
    Jtemp(3,:) = [r s t ones(4,1)]\z(:,k);
    
    trasl(:,k) = Jtemp(:,4);
    J(:,:,k) = Jtemp(:,1:3);
    Jdet(k) = det(J(:,:,k));
    Jcof(:,:,k) = Jdet(k)*inv(J(:,:,k)');
    
end

end

