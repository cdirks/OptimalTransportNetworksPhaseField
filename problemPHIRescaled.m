
function A = problemPHIRescaled ( nodes,elements,sigma,epsilon,phi,alpha,beta )

global BinvDQuad;
global rows;
global cols;

n = length(nodes);
numPhaseFields = size(phi,2);
valsSPhi = zeros(9*size(elements,1),numPhaseFields);
valsMPhi = valsSPhi;

I = eye(2);
D = [-1 -1;I]';

quadPoints = [1 1]'/3;
valCoeffs = [1-quadPoints(1,:)-quadPoints(2,:);quadPoints(1,:);quadPoints(2,:)]';

counter = 1;

for k = 1:size(elements,1)
    
    nodeInds = elements(k, 1:3);    % get global indices
    x = nodes(nodeInds, 1:2);       % local triangle vertices (as matrix rows)
    B = (D * x)';                   % Jacobian of coordinate transform from T_ref to T.
    vol = abs(det(B)/2);            % triangle-volume
    sigmaVals = sigma(k,:);
    sigmaNrmSqr = sigmaVals(:,1).^2+sigmaVals(:,2).^2;
    phiVals = valCoeffs*phi(nodeInds,:);
    phiValsSqrEta = phiVals.^2+alpha*epsilon^2;
    region = find ( phiValsSqrEta == min(phiValsSqrEta), 1, 'first' );
    
    % assemble matrix   
    for i = 1:numPhaseFields 
        vals = vol*beta(i)^2*epsilon/(2*alpha(i))*BinvDQuad{k};
        valsSPhi(counter:counter+8,i) = vals(:);
        vals = [2,1,1;1,2,1;1,1,2] * (vol/12) * ( 2*sigmaNrmSqr/epsilon*(region==i) + alpha(i)/(2*epsilon) );
        valsMPhi(counter:counter+8,i) = vals(:);
    end
        
    counter = counter + 9;
    
end

A = [];
for i = 1:numPhaseFields 
    S = sparse(rows,cols,valsSPhi(:,i),n,n);
    M = sparse(rows,cols,valsMPhi(:,i),n,n);
    A_aux = S+M;
    A = blkdiag(A,A_aux);
end

end