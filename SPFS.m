
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Single phase field solver                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma,phi,lambda] = SPFS(n,epsilonStart,epsilonEnd,maxIter,tol,a,b,example,savename)
disp('+++++ One phase field, no diffuse component +++++');
%-------------------------------------------------------------------------
% Define grid 
x = linspace(0,1,n);
[X,Y] = meshgrid(x,x);
nodes = [X(:) Y(:)];
XAux = X; 
YAux = Y;
elements = delaunay( [XAux(:),YAux(:)] );
numElements = size(elements,1);
numNodes = n^2;
numPhaseFields = length(a);
initialize(nodes,elements);
initialValues = [ones(numNodes*numPhaseFields,1);zeros(numElements*2,1)];
%-------------------------------------------------------------------------
% Define constraint 
feps = getExample(example,n);
SmoothKernel = exp(-(x-.5).^2/(2*(epsilonEnd/2)^2));
SmoothKernel = fftshift(SmoothKernel) / sum(SmoothKernel);
feps = real(ifft2(fft2(feps).*fft2(SmoothKernel'*SmoothKernel)));
feps = transform(feps,nodes,1);
[M,~] = massStiffMatrix(nodes,elements);
[massCentersX,massCentersY] = defineMassCenters(nodes,elements); 
b_aux = M*ones(numNodes,1);
%-------------------------------------------------------------------------
% Equality constraint 
M1 = mixedMassStiffMatrix(nodes,elements,1);
M2 = mixedMassStiffMatrix(nodes,elements,2);
Aeq = [M1' M2'];
beq = -M*feps(:);
%-------------------------------------------------------------------------
% Boundary conditions
boundVec = zeros(numNodes,1);
for k = 1:size(nodes,1)
    if ( nodes(k,1) == 0 || nodes(k,2) == 0 || nodes(k,1) == 1 || nodes(k,2) == 1 ) 
        boundVec(k) = 1; 
    end
end
boundInds = find(boundVec==1);
%-------------------------------------------------------------------------
% Initialization
phi = reshape(initialValues(1:numNodes*numPhaseFields),numNodes,numPhaseFields);
sigma = reshape(initialValues(numNodes*numPhaseFields+1:numNodes*numPhaseFields+numElements*2),numElements,2);
iter = 0;
%-------------------------------------------------------------------------
% ITERATION
figure();
startIter = 20;
endIter = 20;
while ( iter < maxIter+startIter+endIter )
    phiOld = phi;
    sigmaOld = sigma;
    % Set epsilon and eta
    if ( iter > startIter && iter <= maxIter+startIter ) 
        epsilon = epsilonStart-(epsilonStart-epsilonEnd)/(maxIter-1)*(iter-startIter-1);
    elseif ( iter > maxIter+startIter )
        epsilon = epsilonEnd;
    else
        epsilon = epsilonStart;
    end
    % Solve for phi 
    APHI = problemPHI ( nodes,elements,sigma,epsilon,phi,a,b );
    bPHI = b/epsilon*b_aux;
    % Implement boundary conditions
    bPHI = bPHI - APHI*boundVec;
    APHI(:,boundInds) = 0; APHI(boundInds,:) = 0;
    APHI(sub2ind(size(APHI),boundInds,boundInds)) = 1; 
    bPHI(boundInds) = 0;
    phi = APHI\bPHI + boundVec;
    % Solve for lambda 
    ALAMBDA = matrixLAMBDA ( nodes,elements,phi,a,epsilon );
    bLAMBDA = -M*feps;
    lambda = ALAMBDA\bLAMBDA;
    % Set sigma 
    sigma = setSigma(elements,lambda,phi,a,epsilon);
    % Visualize and save
    visualisation(phi,sigma,massCentersX,massCentersY,a,b,n);
    save(savename,'sigma','phi','lambda','epsilon','a','b');
    % Compute error 
    error = max(abs(phiOld(:)-phi(:))) + max(abs(sigmaOld(:)-sigma(:)));
    disp(['Iteration ',num2str(iter),': Error = ',num2str(error)]);
    if ( error < tol && iter < startIter )
        iter = startIter;
    elseif ( error < tol && iter > maxIter+startIter ) 
        break; 
    else
        iter = iter + 1;
    end   
    % Check if constraint is satisfied 
    constraint = Aeq*sigma(:)-beq;
    constraintError = max(abs(constraint(:)));
    disp(['Constraint error = ',num2str(constraintError)]);
end
%-------------------------------------------------------------------------
end