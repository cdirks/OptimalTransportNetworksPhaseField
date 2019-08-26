
function [Mxi,J] = SpecialWeightedMassMatrixAndJacobian(nodes,elements,phi,sigma,alpha0,eta,epsilon,M1,M2)

global valCoeffs;

numElements = size(elements,1);
numNodes = size(nodes,1);
vals = zeros(numElements,1);
valsJ11 = vals;
valsJ12 = vals;
valsJ22 = vals;

I = eye(2);
D = [-1 -1;I]';

for k = 1:numElements
    nodeInds = elements(k,1:3);  
    x = nodes(nodeInds,1:2);        % local triangle vertices (as matrix rows)
    B = (D * x)';                   % Jacobian of coordinate transform from T_ref to T.
    vol = abs(det(B)/2);            % triangle-volume
    phiVals = valCoeffs*phi(nodeInds,:);
    phiValsSqrEta = phiVals.^2+eta.^2;
    gammaEps = min(phiValsSqrEta);
    sigmaNrm = sqrt(sigma(k,1)^2+sigma(k,2)^2);
    % Check cases
    if ( sigmaNrm <= alpha0/(2*gammaEps/epsilon) )
        vals(k) = vol*2*gammaEps/epsilon;
        valsJ11(k) = vals(k);
        valsJ12(k) = 0;
        valsJ22(k) = vals(k);        
    else 
        vals(k) = vol*alpha0/sigmaNrm;
        valsJ11(k) = vals(k) - vol*alpha0*sigma(k,1)^2/(sigmaNrm^3);
        valsJ12(k) = -vol*alpha0*sigma(k,1)*sigma(k,2)/(sigmaNrm^3);
        valsJ22(k) = vals(k) - vol*alpha0*sigma(k,2)^2/(sigmaNrm^3);
    end
end

Mxi = sparse(1:numElements,1:numElements,vals,numElements,numElements);

J11 = sparse(1:numElements,1:numElements,valsJ11,numElements,numElements);
J12 = sparse(1:numElements,1:numElements,valsJ12,numElements,numElements);
J22 = sparse(1:numElements,1:numElements,valsJ22,numElements,numElements);
I1 = sparse(numNodes,numNodes); 

J = [J11,J12,-M1; J12,J22,-M2; M1',M2',I1];

end

