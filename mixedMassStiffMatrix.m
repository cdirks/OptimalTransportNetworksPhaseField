% nodes = n by 2 array whose rows are the mesh nodes
% elements = m by 3 array whose rows are the indices of the nodes of a triangle
% k is coordinate direction for derivative

function S = mixedMassStiffMatrix(nodes,elements,l)
% computes mixed mass and stiffness matrix ( \int_\Omega \phi_i \partial_l \phi_j dx )_{ij}

  n = length(nodes);
  m = size(elements,1);
  rows = zeros(3*size(elements,1),1);
  cols = zeros(3*size(elements,1),1);
  vals = rows;

  I = eye(2);
  D = [-1 -1;I]';

  counter = 1;
  for k = 1:size(elements,1)          % loop through all elements {T}
      nodeInds = elements(k, 1:3);    % get global indices
      x = nodes(nodeInds, 1:2);       % local triangle vertices (as matrix rows)
      A = (D * x)';                   % Jacobian of coordinate transform from T_ref to T.
      AinvD = (D'/A)';                % gradient operator
      vol = abs(det(A)/2);            % triangle-volume

      % assemble matrices
      rows(counter:counter+2) = k;
      cols(counter:counter+2) = nodeInds;
      vals(counter:counter+2) = vol*AinvD(l,:);
      counter = counter + 3;
  end

  S = sparse(rows,cols,vals,m,n);
end
