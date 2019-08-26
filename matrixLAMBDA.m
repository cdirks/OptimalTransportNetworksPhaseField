
function A = matrixLAMBDA ( nodes,elements,phi,alpha,epsilon ) 

global BinvDQuad;
global rows;
global cols;

n = length(nodes);
numElements = size(elements,1);
vals = zeros(4*length(elements),1);

quadPoints = [1 1]'/3;
valCoeffs = [1-quadPoints(1,:)-quadPoints(2,:);quadPoints(1,:);quadPoints(2,:)]';

I = eye(2);
D = [-1 -1;I]';

counter = 1;

for k = 1:numElements
    
    nodeInds = elements(k,1:3);     % get global indices
    x = nodes(nodeInds,1:2);        % local triangle vertices (as matrix rows)
    B = (D * x)';                   % Jacobian of coordinate transform from T_ref to T.
    vol = abs(det(B)/2);            % triangle-volume
    phiVals = valCoeffs*phi(nodeInds,:);   
    phiValsSqrEta = phiVals.^2 + alpha*epsilon^2;
    gammaEps = min(phiValsSqrEta);
    
    % assemble matrix
    vals(counter:counter+8) = vol*epsilon/(2*gammaEps)*BinvDQuad{k};
    
    counter = counter + 9;
    
end

A = sparse(rows,cols,vals,n,n);

end