function [x] = solve_linear_system_with_LU(M, rhs)

%preprocess
    [M_data.L, M_data.U, M_data.p, M_data.q] = lu(M, 'vector');
%solve
	x(M_data.q, :) = M_data.U\(M_data.L\(rhs(M_data.p, :)));
end
