function [nod2, wei2, nod3, wei3, node_maps, node_maps_inv] = quadrature(N, Nfaces)
%[nod2, wei2, nod3, wei3, node_maps, node_maps_inv] = quadrature(N, Nfaces).
%Computes quadrature nodes and weights over the simplex tetrahedron
%and triangle. 3D from http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf and 2D from
%Riviere's book; commented ones from Quarteroni's book.
%Computes also the maps from the 2D refence triangle to the faces of
%the 3D reference tetrahedron.
%Ones in the quadrature nodes are needed by the translation .

if N == 1
    % Degree of exactness 2
    nod3 = [0.585410196624969 0.138196601125011 0.138196601125011 0.138196601125011;
            0.138196601125011 0.585410196624969 0.138196601125011 0.138196601125011;
            0.138196601125011 0.138196601125011 0.585410196624969 0.138196601125011;
            1                 1                 1                 1];
    wei3 = 1/6*[0.25 0.25 0.25 0.25];
    % Degree of exacness 2
    nod2 = [2/3 1/6 1/6;
            1/6 1/6 2/3;
             1   1   1;
             1   1   1];
    wei2 = 0.5*[1/3 1/3 1/3];
    
elseif N == 2
    % Degree of exactness 2
    nod3 = [0.585410196624969 0.138196601125011 0.138196601125011 0.138196601125011;
            0.138196601125011 0.585410196624969 0.138196601125011 0.138196601125011;
            0.138196601125011 0.138196601125011 0.585410196624969 0.138196601125011;
            1                 1                 1                 1];
    wei3 = 1/6*[0.25 0.25 0.25 0.25];
    % Degree of exactness 4
    nod2 = [0.108103018168 0.445948490915 0.445948490915 0.816847572980 0.091576213509 0.091576213509;
            0.445948490915 0.445948490915 0.108103018168 0.091576213509 0.091576213509 0.816847572980;
            1              1              1              1              1              1;
            1              1              1              1              1              1];
    wei2 = [0.1116907948390 0.1116907948390 0.1116907948390 0.0549758718276 0.0549758718276 0.0549758718276];
    
elseif N >= 3
    % Degree of exactness 5
    nod3 = [1/4 0 1/3 1/3 1/3 8/11 1/11 1/11 1/11 0.4334498464263357 0.0665501535736643 0.0665501535736643 0.0665501535736643 0.4334498464263357 0.4334498464263357;
            1/4 1/3 1/3 1/3 0 1/11 1/11 1/11 8/11 0.0665501535736643 0.4334498464263357 0.0665501535736643 0.4334498464263357 0.0665501535736643 0.4334498464263357;
            1/4 1/3 1/3 0 1/3 1/11 1/11 8/11 1/11 0.0665501535736643 0.0665501535736643 0.4334498464263357 0.4334498464263357 0.4334498464263357 0.0665501535736643;
            1   1   1  1  1   1    1    1    1    1                  1                  1                  1                  1                  1];
    wei3 = [0.1817020685825351 0.0361607142857143 0.0361607142857143 0.0361607142857143 0.0361607142857143 ...
            0.0698714945161738 0.0698714945161738 0.0698714945161738 0.0698714945161738 0.0656948493683187 ...
            0.0656948493683187 0.0656948493683187 0.0656948493683187 0.0656948493683187 0.0656948493683187]/6;
    % Degree of exactness 7
    nod2 = [1/3 0.479308067841 0.260345966079 0.260345966079 0.869739794195 0.065130102902 0.065130102902 0.048690315425 0.312865496004 0.048690315425 0.638444188569 0.638444188569 0.312865496004;
        1/3 0.260345966079 0.260345966079 0.479308067841 0.065130102902 0.065130102902 0.869739794195 0.312865496004 0.048690315425 0.638444188569 0.048690315425 0.312865496004 0.638444188569;
        1   1              1              1              1              1              1              1              1              1              1              1              1;
        1   1              1              1              1              1              1              1              1              1              1              1              1];    
    wei2 = [-0.074785022233 0.087807628716 0.087807628716 0.087807628716 0.026673617804 0.026673617804 0.026673617804 0.038556880445 0.038556880445 0.038556880445 0.038556880445 0.038556880445 0.038556880445];    
end
    
% Degree of exactness 5  
% nod3 = [0.25 0.5 1/6 1/6 1/6;
%         0.25 1/6 0.5 1/6 1/6;
%         0.25 1/6 1/6 0.5 1/6;
%          1    1   1   1   1];
% wei3 = 1/6*[-4/5 9/20 9/20 9/20 9/20];

% Degree of exactness 3
% nod2 = [1/3 0.6 0.2 0.2;
%         1/3 0.2 0.6 0.2;
%          1   1   1   1;
%          1   1   1   1]; % ones are needed by the transaltion 
% wei2 = 0.5*[-9/16 25/48 25/48 25/48];

% Maps from the 2D refence triangle to the faces of the 3D reference
% tetrahedron and back.
node_maps = zeros(4,4,Nfaces);
node_maps_inv = zeros(2,3,Nfaces);
node_maps(:,:,1) = [0 1 0 0;
                    1 0 0 0;
                    0 0 0 0;
                    0 0 0 1];
             
node_maps(:,:,2) = [1 0 0 0;
                    0 0 0 0;
                    0 1 0 0;
                    0 0 0 1];
             
node_maps(:,:,3) = [0 0 0 0;
                    0 1 0 0;
                    1 0 0 0;
                    0 0 0 1];
             
node_maps(:,:,4) = [1  0 0 0;
                    0  1 0 0;
                   -1 -1 1 0;
                    0  0 0 1];
               
node_maps_inv(:,:,1) = [0 1 0;
                        1 0 0];
                    
node_maps_inv(:,:,2) = [1 0 0;
                        0 0 1];

node_maps_inv(:,:,3) = [0 0 1;
                        0 1 0];
           
node_maps_inv(:,:,4) = [1 0 0;
                        0 1 0];
              
end