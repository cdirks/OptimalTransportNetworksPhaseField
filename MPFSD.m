
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Multiple phase field solver with diffuse component          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma,phi,lambda] = MPFSD(n,epsilonStart,epsilonEnd,maxIter,tol,a0,a,b,example,savename,initSol,savenameInit)
numPhaseFields = length(a);
disp(['+++++ ',num2str(numPhaseFields),' phase fields with diffuse component +++++']);
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
% Apply mass factor 
massFactor = 1;
a = a/massFactor;
feps = feps*massFactor;
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
boundIndsMult = zeros(numPhaseFields*length(boundInds),1);
for k = 1:numPhaseFields
    boundIndsMult((k-1)*length(boundInds)+1:k*length(boundInds)) = boundInds+(k-1)*numNodes; 
end
boundVecMult = repmat(boundVec,numPhaseFields,1);
%-------------------------------------------------------------------------
% Define phi on elements 
nodeToElement = sparse([1:numElements 1:numElements 1:numElements],[elements(:,1) elements(:,2) elements(:,3)],[1/3*ones(numElements,1) 1/3*ones(numElements,1) 1/3*ones(numElements,1)],numElements,numNodes);
aux = 1./sum(nodeToElement);
elementToNode = spdiags(aux',0,numNodes,numNodes)*(nodeToElement');
%-------------------------------------------------------------------------
% Initialization
phi = reshape(initialValues(1:numNodes*numPhaseFields),numNodes,numPhaseFields);
sigma = reshape(initialValues(numNodes*numPhaseFields+1:numNodes*numPhaseFields+numElements*2),numElements,2);
%-------------------------------------------------------------------------
% Compute initial sigma
if ( ~initSol ) 
    startIter = 0;
    endIter = 0;
    iter = 0;
    while ( iter < startIter+maxIter+endIter )       
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
        % Solve problem in phi
        APHI = problemPHIRescaled ( nodes,elements,sigma,epsilon,phi(:,1),a(1),b(1) );
        bPHI = a(1)/(2*epsilon)*b_aux; 
        % Implement boundary conditions
        bPHI = bPHI - APHI*boundVec;
        APHI(:,boundInds) = 0; APHI(boundInds,:) = 0;
        APHI(sub2ind(size(APHI),boundInds,boundInds)) = 1; 
        bPHI(boundInds) = 0;
        phi(:,1) = APHI\bPHI + boundVec;
        % Solve problem in lambda 
        ALAMBDA = matrixLAMBDA ( nodes,elements,phi(:,1),a(1),epsilon );
        bLAMBDA = -M*feps;
        lambda = ALAMBDA\bLAMBDA;
        % Set sigma 
        sigma = setSigma(elements,lambda,phi(:,1),a(1),epsilon);
        % Visualize
        visualisation(phi(:,1),sigma,massCentersX,massCentersY,a(1),b(1),n);
        % Compute error 
        error = max(abs(phiOld(:,1)-phi(:,1))) + max(abs(sigmaOld(:)-sigma(:)));
        save(savenameInit,'sigma','lambda');
        disp(['Initial iteration ',num2str(iter),': Error = ',num2str(error)]);
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
else
    load(savenameInit,'sigma'); 
    load(savenameInit,'lambda');
end
%-------------------------------------------------------------------------
% Set initial phase fields   
phi = ones(numNodes,numPhaseFields+1);
phiTest = phi;
% Set epsilon
epsilon = epsilonEnd;
% Solve for every phi independently 
for i = 1:numPhaseFields
    APHI = problemPHIRescaled ( nodes,elements,sigma,epsilon,phi(:,i),a(i),b(i) );
    bPHI = kron(a(i)/(2*epsilon),b_aux); 
    phiTest(:,i+1) = APHI\bPHI;
end
phiTest = reshape(phiTest,numNodes,numPhaseFields+1);
% Compute mass from sigma
width = max(8/(n-1),2*epsilon); 
hx = 1/(n-1)^2;
sigmaNode = elementToNode*sigma;
sigmaNodeNrm = sqrt(sigmaNode(:,1).^2+sigmaNode(:,2).^2);
sigmaNodeNrm = reshape(sigmaNodeNrm,n,n);
sigmaNodeNrm = sigmaNodeNrm';
sigmaNodeNrm = sigmaNodeNrm(:);
sigmaConv = zeros(n,n);
for i = 1:n
    for j = 1:n
        xAux = [x(i)-X(:) x(j)-Y(:)];
        GaussKernel = double(sqrt(xAux(:,1).^2+xAux(:,2).^2)<=width);
        integral = hx/(2*width) * sigmaNodeNrm'*GaussKernel;
        sigmaConv(i,j) = integral;
    end    
end 
% Compute costs depending on mass
nodeEnergy = zeros(numNodes,numPhaseFields+1);
nodeEnergy(:,1) = a0*sigmaConv(:);
for i = 1:numPhaseFields 
    nodeEnergy(:,i+1) = b(i)+sigmaConv(:)*a(i); 
end
% Compute locally optimal phase fields
[~,minIndex] = min(nodeEnergy,[],2);
linearIndex = [(1:numNodes)' minIndex];
linearIndex = linearIndex(:,1)+(linearIndex(:,2)-1)*numNodes;
phi(linearIndex) = phiTest(linearIndex);
phi = phi(:,2:end);
%-------------------------------------------------------------------------
% ITERATION
figure();
startIter = 20;
endIter = 20;
iter = 0;
error = 1;
I1 = sparse(numElements,numElements);
I2 = sparse(numNodes,numNodes);
I3 = sparse(2*numElements,1);
while ( iter < maxIter+startIter+endIter && error > tol )
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
    eta = sqrt(a)*epsilon; 
    % Update phi (explicit)
    APHI = problemPHIDiffuse ( nodes,elements,sigma,epsilon,phi,a0,a,b );
    bPHI = kron(a(:)./(2*epsilon),b_aux); 
    % Implement boundary conditions
    bPHI = bPHI - APHI*boundVecMult; % TODO: Hier weiter
    APHI(:,boundIndsMult) = 0; APHI(boundIndsMult,:) = 0;
    APHI(sub2ind(size(APHI),boundIndsMult,boundIndsMult)) = 1; 
    bPHI(boundIndsMult) = 0;
    phi = APHI\bPHI + boundVecMult;
    phi = reshape(phi,numNodes,numPhaseFields);
    % Update lambda and sigma (Newton step) 
    [Mxi,J] = SpecialWeightedMassMatrixAndJacobian(nodes,elements,phi,sigma,a0,eta,epsilon,M1,M2);        
    R = [Mxi,I1,-M1; I1,Mxi,-M2; M1',M2',I2] * [sigma(:);lambda(:)] + [I3;M*feps(:)];
    sol = gmres(J,R); 
    sigma(:,1) = sigma(:,1) - sol(1:numElements);
    sigma(:,2) = sigma(:,2) - sol(numElements+1:2*numElements);
    lambda = lambda - sol(2*numElements+1:end);
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