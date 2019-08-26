%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Visualisation tool                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualisation(phi,sigma,massCentersX,massCentersY,a,b,n)

% phi
numPhaseFields = size(phi,2);
for i = 1:numPhaseFields 
    subplot(1,numPhaseFields+1,i);
    imagesc(reshape(phi(:,i),n,n));
    title(['a=',num2str(a(i)),', b=',num2str(b(i))]);
    axis image; 
    drawnow;
end

% sigma
subplot(1,numPhaseFields+1,numPhaseFields+1);
quiver(massCentersX,massCentersY,sigma(:,1),sigma(:,2));
title('Flux');
axis image;
drawnow;

end





