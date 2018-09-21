
function Hamiltonian = kraken_1d(Nx, dx, V_1D, V_int)

% http://www.mathworks.com/help/matlab/ref/sparse.html
% see S = sparse(i,j,v)
% these vectors have 3*Nx-2 elements because we have Nx elements in the
% diagonal and Nx-1 in each of the sub-/super- diagonals

sparse_elements = 3*Nx - 2;

sparse_i = zeros(1, sparse_elements); % sparse matrix i-indeces
sparse_j = zeros(1, sparse_elements); % sparse matrix j-indeces
sparse_v = zeros(1, sparse_elements); % sparse matrix values at positions (i,j) 

% Nx DIAGONAL ELEMENTS
for i=1:Nx
    
    index = i;
     
    sparse_i(index) = i;
    sparse_j(index) = i;
    sparse_v(index) = V_1D(i) + 1./dx^2;
    
end

% Nx-1 ABOVE THE DIAGONAL ELEMENTS
for i = 1 : Nx-1
    
    index = Nx + i;
     
    sparse_i(index) = i+1;
    sparse_j(index) = i;
    sparse_v(index) = - 0.5 / dx^2;

end

% Nx-1 BELOW THE DIAGONAL ELEMENTS
for i = 1 : Nx-1
    
    index = Nx + Nx-1 + i;

    sparse_i(index) = i;
    sparse_j(index) = i+1;
    sparse_v(index) = - 0.5 / dx^2;
    
end


%% RELEASE THE KRAKEN!
Hamiltonian = sparse(sparse_i, sparse_j, sparse_v);
