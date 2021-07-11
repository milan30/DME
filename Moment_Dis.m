function mom = Moment_Dis_up(m,N,P,D,pos)
if nargin<5
    pos = N/2;
end
z=repmat(1:5,5,1);
temp_mat = tril(z,-1)-tril(z,-2);
B = temp_mat';
B_square_by_2 = 1/2.*B^2;
B_cube_by_6 = 1/6.*B^3;
B_fourth_by_24 = 1/24.*B^4;
DP = D*P;
D_sq_P = D^2*P;
D_cube_P = D^3*P;
D_fourth_P = D^4*P;
der_vec = zeros(N,5);
der_vec(:,1)=ones(N,1);
for j=1:(m+1)
    der_vec = P*der_vec - DP*der_vec* B + D_sq_P*der_vec*B_square_by_2 - D_cube_P*der_vec*B_cube_by_6...
    + D_fourth_P*der_vec*B_fourth_by_24; 
end

%%finding specific derivative at pos-th position
vec_N_by_2 = zeros(1,N);
vec_N_by_2(pos) = 1; 
mom = vec_N_by_2* der_vec;
mom = mom(2:5);
end
