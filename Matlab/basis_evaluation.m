function [phi, dphi, phi_bordo, grad_bordo, phi_func] = basis_evaluation(nod3, nod2, node_maps, N, Np)
%[phi, dphi, phi_bordo, grad_bordo, phi_func] = basis_evaluation(N, nod3, nod2, node_maps, N, Np)
%Gives the evaluation of the P1 or P2 basis functions and their gradient
%at the 3D quadrature nodes and 2D face quadrature nodes.

if N == 1
    phi_func = {@(p) (-p(1,:)-p(2,:)-p(3,:)+1);
                @(p) (p(1,:));
                @(p) (p(2,:));
                @(p) (p(3,:))};

    grad_phi_func = {@(p) [-1;-1;-1];
                     @(p) [1;0;0];
                     @(p) [0;1;0];
                     @(p) [0;0;1]};

elseif N == 2
    phi_func = {@(p) 2*(p(1,:)+p(2,:)+p(3,:)).^2-3*(p(1,:)+p(2,:)+p(3,:))+1;
                @(p) 2*p(1,:).^2-p(1,:);
                @(p) 2*p(2,:).^2-p(2,:);
                @(p) 2*p(3,:).^2-p(3,:);
                @(p) -4*p(1,:).*(p(1,:)+p(2,:)+p(3,:)-1);
                @(p) -4*p(2,:).*(p(1,:)+p(2,:)+p(3,:)-1);
                @(p) -4*p(3,:).*(p(1,:)+p(2,:)+p(3,:)-1);
                @(p) 4*p(1,:).*p(2,:);
                @(p) 4*p(1,:).*p(3,:);
                @(p) 4*p(2,:).*p(3,:)};
            
    grad_phi_func = {@(p) [4*(p(1,:)+p(2,:)+p(3,:))-3; 4*(p(1,:)+p(2,:)+p(3,:))-3; 4*(p(1,:)+p(2,:)+p(3,:))-3];
                     @(p) [4*p(1,:)-1; 0; 0];
                     @(p) [0; 4*p(2,:)-1; 0];
                     @(p) [0; 0; 4*p(3,:)-1];
                     @(p) [-4*(2*p(1,:)+p(2,:)+p(3,:)-1); -4*p(1,:); -4*p(1,:)];
                     @(p) [-4*p(2,:); -4*(p(1,:)+2*p(2,:)+p(3,:)-1); -4*p(2,:)];
                     @(p) [-4*p(3,:); -4*p(3,:); -4*(p(1,:)+p(2,:)+2*p(3,:)-1)];
                     @(p) [4*p(2,:); 4*p(1,:); 0];
                     @(p) [4*p(3,:); 0; 4*p(1,:)];
                     @(p) [0; 4*p(3,:); 4*p(2,:)]};
            
else
    disp('Wrong degree');
end

nq3 = size(nod3, 2);
nq2 = size(nod2, 2);

%three derivatives, nq3 quadrature points and Np functions
dphi = zeros(3,nq3,Np);
phi = zeros(nq3,Np);
    
%evaluation on 3 derivatives, nq2 quadrature boundary points, 4 faces, Np functions.
grad_bordo = zeros(3,nq2,4,Np);
phi_bordo = zeros(nq2,4,Np);

for f = 1:Np
    phi(:,f) = phi_func{f}(nod3)';
    for q = 1:nq3
        dphi(:,q,f) = grad_phi_func{f}(nod3(:,q));
    end
    
    
    %dphi(:,:,f) = repmat(grad_phi_func(:,f), [1 nq3]);
    for e = 1:4
        for q = 1:nq2
            grad_bordo(:,q,e,f) = grad_phi_func{f}(node_maps(:,:,e)*nod2(:,q));
        end
        phi_bordo(:,e,f) = phi_func{f}(node_maps(:,:,e)*nod2);
    end
end

end
