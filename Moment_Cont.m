function mom= Moment_Cont(t, N, G, D,pos)
if nargin<5
    pos = N/2;
end
% S is a diagonal matrix of order N
eye_kron_G = kron(eye(5),G);
% constructing B transpose
B_con=repmat(1:5,5,1);
B_transpose = tril(B_con,-1)-tril(B_con,-2);
B_kron_S = kron(B_transpose, D);
main_mat = eye_kron_G - B_kron_S;
initial_con = kron(eye(5,1),ones(N,1));
%using Pade approximation
sol_ode = expm(t*main_mat)*initial_con;
%using eigen decomposition
%sol_ode = funm(t*main_mat,@exp)*initial_con;
%using Krylov approximation
%sol_ode = expv(t,main_mat,initial_con);
%using Markov algorithm
%sol_ode = mexpv(t,main_mat,initial_con);
%using Chebyshev algorithm
%sol_ode = chbv(main_mat,initial_con);
%finding specific derivative at pos-th position
vec_N_by_2 = zeros(1,N);
vec_N_by_2(pos) = 1; 
mom = kron(eye(5),vec_N_by_2)*sol_ode;
mom = mom(2:5);


%time consuming approach
% sol_ode = reshape(sol_ode,N,5); 
% mom1 = sol_ode(:,2:5);
% mom = mom1(N/2,:);
end
