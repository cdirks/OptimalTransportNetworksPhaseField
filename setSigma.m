
function sigma = setSigma(elements,lambda,phi,alpha,epsilon)

global BinvD;
global valCoeffs;

numElements = size(elements,1);
sigma = zeros(numElements,2);

for k = 1:size(elements,1)
    
    nodeInds = elements(k,1:3);                        % get global node indices
    lambdaGrad = BinvD{k} * lambda(nodeInds,:);        % gradient of the lambda on element
    
    phiVals = valCoeffs*phi(nodeInds,:);
    phiValsSqrEta = phiVals.^2+alpha*epsilon^2;
    gammaEps = min(phiValsSqrEta);
    
    sigma(k,1) = epsilon*lambdaGrad(1,:)./(2*gammaEps);
    sigma(k,2) = epsilon*lambdaGrad(2,:)./(2*gammaEps);

end

end