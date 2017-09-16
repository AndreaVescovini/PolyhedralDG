function [P] = LegendrePder(x, N, int)
% function [P] = LegendreP(x, N, int)
% Purpose: Evaluate the first derivative of Legendre Polynomial in int
% at points x for order N and returns P[1:length(xp))]

% Turn points into row if needed.
xp = x;
dims = size(xp);
if (dims(2) == 1)
    xp = xp';
end

PL = zeros(N+1,length(xp));
P = zeros(1,length(xp));

hb = (int(2)-int(1))/2;
mb = (int(1)+int(2))/2;
xp = (xp - mb)/hb;

% Initial values L_0(x) and L_1(x)
PL(1,:) = 1.0;
if (N==0)
    P = (P')*(1/hb)^1.5;
    return;
end
PL(2,:) = xp;

% Forward recurrence using the symmetry of the recurrence.
for i = 2:N-1
  PL(i+1,:) = -(i-1)/i*PL(i-1,:) + (2*i-1)/i*xp.*PL(i,:);
end

for k = 0:floor((N-1)/2)
    P = P + PL(N-2*k,:)*(2*N-4*k-1);
end

P = (P')*(1/hb)^1.5;
end