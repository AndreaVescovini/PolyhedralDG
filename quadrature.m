function [nod2, wei2, nod3, wei3] = quadrature()
%quadrature. Computes quadrature nodes and weights over the simplex tetrahedron,
%degree of exactness 3 and quadrature nodes and weights over the simplex triangle,
%degree of exactness 3 (Modellistica numerica per problemi differenziali,
%Quarteroni, pag 180).
   
nod3 = [-0.5   0  -2/3 -2/3 -2/3;
        -0.5 -2/3   0  -2/3 -2/3;
        -0.5 -2/3 -2/3   0  -2/3];
   
wei3 = 8/3*[-4/5; 9/20; 9/20; 9/20; 9/20];

nod2 = [-1/3 1/5 -3/5 1/5;
        -1/3 -3/5 1/5 1/5];
     
wei2 = 2*[-9/16; 25/48; 25/48; 25/48];
end