
function Hamiltonian = kraken_3d(Nx, dx, V_1D, V_int)

% http://www.mathworks.com/help/matlab/ref/sparse.html
% see S = sparse(i,j,v)
% these vectors have 3*Nx-2 elements because we have Nx elements in the
% diagonal and Nx-1 in each of the sub-/super- diagonals

sparse_elements = 7*Nx^3 - 6*Nx^2;

sparse_i = zeros(1, sparse_elements); % sparse matrix i-indeces
sparse_j = zeros(1, sparse_elements); % sparse matrix j-indeces
sparse_v = zeros(1, sparse_elements); % sparse matrix values at positions (i,j) 

s_index = 0;

% Nx*Nx*Nx DIAGONAL ELEMENTS
for i = 1 : Nx
    for j = 1 : Nx
        for k = 1 : Nx
    
            k_index = (i-1)*Nx*Nx + (j-1)*Nx + k;
            s_index = s_index + 1;

            sparse_i(s_index) = k_index;
            sparse_j(s_index) = k_index;
            sparse_v(s_index) = V_1D(i) + V_1D(j) + V_1D(k) + V_int(i,j) + V_int(i,k) + V_int(j,k) + 1./dx^2 + 1./dx^2 + 1./dx^2;% + g*kroneckerDelta(i,j);
            
        end    
    end
end

% Nx*Nx*(Nx-1) ABOVE THE DIAGONAL ELEMENTS
for i = 1 : Nx
    for j = 1 : Nx
        for k = 1 : Nx-1
    
            k_index = (i-1)*Nx*Nx + (j-1)*Nx + k;
            s_index = s_index + 1;

            sparse_i(s_index) = k_index + 1;
            sparse_j(s_index) = k_index;
            sparse_v(s_index) = - 0.5 / dx^2;
            
        end
    end
end

% Nx*Nx*(Nx-1) BELOW THE DIAGONAL ELEMENTS
for i = 1 : Nx
    for j = 1 : Nx
        for k = 1 : Nx-1
    
            k_index = (i-1)*Nx*Nx + (j-1)*Nx + k;
            s_index = s_index + 1;

            sparse_i(s_index) = k_index;
            sparse_j(s_index) = k_index + 1;
            sparse_v(s_index) = - 0.5 / dx^2;

        end
    end
end

% Nx*(Nx-1)*Nx SECOND ABOVE THE DIAGONAL ELEMENTS
for i = 1 : Nx
    for j = 1 : Nx-1
        for k = 1 : Nx
    
            k_index = (i-1)*Nx*Nx + (j-1)*Nx + k;
            s_index = s_index + 1;     
            
            sparse_i(s_index) = k_index + Nx;
            sparse_j(s_index) = k_index;
            sparse_v(s_index) = - 0.5 / dx^2;

        end
    end
end

% Nx*(Nx-1)*Nx BELOW THE DIAGONAL ELEMENTS
for i = 1 : Nx
    for j = 1 : Nx-1
        for k = 1 : Nx
    
            k_index = (i-1)*Nx*Nx + (j-1)*Nx + k;
            s_index = s_index + 1;     
            
            sparse_i(s_index) = k_index;
            sparse_j(s_index) = k_index + Nx;
            sparse_v(s_index) = - 0.5 / dx^2;

        end
    end
end

% (Nx-1)*Nx*Nx BELOW BELOW THE DIAGONAL ELEMENTS
for i = 1 : Nx-1
    for j = 1 : Nx
        for k = 1 : Nx
    
            k_index = (i-1)*Nx*Nx + (j-1)*Nx + k;
            s_index = s_index + 1;     
            
            sparse_i(s_index) = k_index + Nx*Nx;
            sparse_j(s_index) = k_index;
            sparse_v(s_index) = - 0.5 / dx^2;

        end
    end
end

% (Nx-1)*Nx*Nx ABOVE ABOVE THE DIAGONAL ELEMENTS
for i = 1 : Nx-1
    for j = 1 : Nx
        for k = 1 : Nx
    
            k_index = (i-1)*Nx*Nx + (j-1)*Nx + k;
            s_index = s_index + 1;     
            
            sparse_i(s_index) = k_index;
            sparse_j(s_index) = k_index + Nx*Nx;
            sparse_v(s_index) = - 0.5 / dx^2;

        end
    end
end


%% RELEASE THE KRAKEN!
Hamiltonian = sparse(sparse_i, sparse_j, sparse_v);
