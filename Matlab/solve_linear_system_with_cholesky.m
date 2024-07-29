function [x] = solve_linear_system_with_cholesky(M, rhs)

if ~ishermitian(M)
    x = zeros(size(rhs));
    error("Error: the matrix is not symmetric/hermitian");
    return;
end

%preprocess
    [M_data.R, k, M_data.p] = chol(M, 'vector');
    if(k ~= 0)
        x = zeros(size(rhs));
        error("Error: the matrix is not positive definite");
        return;
    end
    M_data.Rt = M_data.R';
%solve   
    x(M_data.p, :) = (M_data.R\(M_data.Rt\(rhs(M_data.p, :))));

end
