function [D,Dt] = defDDt

D = @(U) ForwardD(U);
Dt = @(X,Y) Dive(X,Y);

function [Dux,Duy] = ForwardD(U)
% [ux,uy] = D u

Dux = [diff(U,1,2), U(:,1) - U(:,end)];
Duy = [diff(U,1,1); U(1,:) - U(end,:)];

function DtXY = Dive(X,Y)
% DtXY = D_1' X + D_2' Y

DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];
DtXY = DtXY(:);