function [P] = LegendreP(x, N, int)
% function [P] = LegendreP(x, N, int)
% Purpose: Evaluate Legendre Polynomial in int at points x for order N and returns P[1:length(xp))]

% Turn points into row if needed.
xp = x;
dims = size(xp);
if (dims(2) == 1)
    xp = xp';
end

PL = zeros(N+1,length(xp));

hb = (int(2)-int(1))/2;
mb = (int(2)+int(1))/2;
xp = (xp - mb)/hb;

% Initial values L_0(x) and L_1(x)
PL(1,:) = 1.0;
if (N==0)
    P=(PL')*sqrt(1/hb);
    return;
end
PL(2,:) = xp;
if (N==1)
    P=(PL(N+1,:)')*sqrt(1/hb);
    return;
end

% Forward recurrence using the symmetry of the recurrence.
for i = 2:N
  PL(i+1,:) = -(i-1)/i*PL(i-1,:) + (2*i-1)/i*xp.*PL(i,:);
end

P = (PL(N+1,:)')*sqrt(1/hb);
end