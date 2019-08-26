% Transform image from FE representation to matrix and vice versa 
function B = transform(A,nodes,dir)
 
n = sqrt(size(nodes,1));
nodesPlot = round(nodes*(n-1)+1);

% FE to matrix
if ( dir == 0 ) 
    B = zeros(n,n);
    for i = 1:n^2
        B(nodesPlot(i,1),nodesPlot(i,2)) = A(i);
    end
    B = reshape(B,n,n);
    B = flip(B);
    B = B';
% Matrix to FE
elseif ( dir == 1 ) 
    Ac = flip(A');
    B = zeros(n^2,1);
    for i = 1:n^2
        B(i) = Ac(nodesPlot(i,1),nodesPlot(i,2)); 
    end
end

end