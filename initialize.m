
function initialize ( nodes,elements ) 

global BinvDQuad;
global BinvD;
global rows;
global cols;
global quadPoints;
global valCoeffs;

I = eye(2);
D = [-1 -1;I]';
BinvDQuad = cell(size(elements,1),1);
BinvD = cell(size(elements,1),1);
rows = zeros(9*length(elements),1);

cols = rows;
counter = 1;

quadPoints = [1 1]'/3;
valCoeffs = [1-quadPoints(1,:)-quadPoints(2,:);quadPoints(1,:);quadPoints(2,:)]';

for k = 1:size(elements,1)
    nodeInds = elements(k, 1:3);    % get global indices
    x = nodes(nodeInds, 1:2);       % local triangle vertices (as matrix rows)
    B = (D * x)';                   % Jacobian of coordinate transform from T_ref to T.
    BinvDTmp = (D'/B)';             % gradient operator
    BinvDQuad{k} = BinvDTmp'*BinvDTmp;
    BinvD{k} = BinvDTmp;
    rows(counter:counter+8) = repmat(nodeInds(:),3,1);
    c = repmat(nodeInds,3,1);
    cols(counter:counter+8) = c(:);
    counter = counter + 9;
end

end