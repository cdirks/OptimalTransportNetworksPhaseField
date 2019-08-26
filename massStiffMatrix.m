% nodes = n by 2 array whose rows are the mesh nodes
% elements = m by 3 array whose rows are the indices of the nodes of a triangle

function [M,L] = massStiffMatrix(nodes,elements)
% computes standard mass matrix M and stiffness matrix L using linear finite elements

  global rows;
  global cols;
  global BinvDQuad;

  n = length(nodes);
  valsM = rows;
  valsL = rows;

  I = eye(2);
  D = [-1 -1;I]';

  counter = 1;
  for k = 1:size(elements,1)          % loop through all elements {T}
      nodeInds = elements(k, 1:3);    % get global indices
      x = nodes(nodeInds, 1:2);       % local triangle vertices (as matrix rows)
      A = (D * x)';                   % Jacobian of coordinate transform from T_ref to T.
      vol = abs(det(A)/2);            % triangle-volume

      % assemble matrices
      valsM(counter:counter+8) = [2,1,1;1,2,1;1,1,2] * (vol/12);
      valsL(counter:counter+8) = vol*BinvDQuad{k};
      counter = counter + 9;
  end

  M = sparse(rows,cols,valsM,n,n);
  L = sparse(rows,cols,valsL,n,n);
end