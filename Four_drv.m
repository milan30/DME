function deri = Four_drv(t, N, A, S)
[P,D] = eig(A);
L =   P\(S * P);
D = diag(D);
exp_D = exp(t*D);
gam1 = zeros(N);
for i=1:N
    for j=1:N
        if D(i) ~= D(j)
            gam1(i,j)=(exp_D(i)-exp_D(j))/(D(i)-D(j));
        else
            gam1(i,j) = t*exp_D(i);
        end
    end
end
gam_tilde_1 = L .* gam1; %theorem 2.9

% portion for second derivative
%gam_3d = zeros(N,N,N);
drv = zeros(N);
gam3d = zeros(N,N,N);

for i=1:N
    temp_mat = zeros(N);
    for j = 1:N
        temp1=0;
        for k=1:N
            if D(k)~=D(j)
                vec_i = (gam1(i,j)-gam1(i,k))/(D(j)-D(k));
            else
                if D(i)~= D(j)
                    vec_i = (gam1(j,j)-gam1(i,j))/(D(j)-D(i));
                else
                    vec_i = t/2* gam1(i,i);
                end
            end

            temp1 = temp1 + L(i,k)*vec_i*L(k,j);
            temp_mat(k,j) = vec_i;
            
        end
        drv(i,j)= temp1;
    end
    gam3d(:,:,i)=temp_mat;
end

dev3=zeros(N);
gam4d = zeros(N,N,N,N);
for i = 1:N
    for j = 1:N
        for l = 1:N
            for k = 1:N
                if (D(l) == D(j))
                    if(D(k) == D(j))
                        if (D(i) == D(j))
                            temp = t^3*exp(t*D(j))/6.0;
                        else
                            temp = (gam3d(i,j,j)-gam3d(j,j,j))/(D(i)-D(j));
                        end
                    else% ! k neq j
                        temp = (gam3d(j,j,i)-gam3d(k,j,i))/(D(j)-D(k));
                    end% IF ! k
                else %! l neq j
                    temp = (gam3d(j,k,i)-gam3d(l,k,i))/(D(j)-D(l));
                end%IF ! l = j
                dev3(i,j) = dev3(i,j) + L(i,k)*L(k,l)*L(l,j)*temp;
                gam4d(k,j,i,l)=temp;
            end
        end
    end
end
        
dev4=zeros(N);
for i = 1:N
    for j = 1:N
        for l = 1:N
            for k = 1: N
                for m = 1:N
                    if D(l) == D(j)
                        if D(k) == D(j)
                            if D(m) == D(j)
                                if D(i) == D(j) % (j, l, k, m, i) = (j, j, j, j, j)
                                    temp1 = t^4*exp(t*D(j))/24;
                                else 	% (j, l, k, m, i) = (j, j, j, j, i)	= (i, j, j, j, j)
                                    temp1 = (gam4d(i,j,j,j)-gam4d(j,j,j,j))/(D(i)-D(j));
                                end %IF ! i
                            else % m neq j, i.e., (j, l, k, m, i) = (j, j, j, m, i) = (m, j, j, j, i)
                                temp1 = (gam4d(m,j,j,i)-gam4d(j,j,j,i))/(D(m)-D(j));
                            end %IF ! m
                        else  % k neq j i.e., (j, l, k, m, i) = (j, j, k, m, i) = (j, k, j, m, i)
                            temp1 = (gam4d(j,j,m,i)-gam4d(k,j,m,i))/(D(j)-D(k));
                        end %IF ! k
                    else %! l neq j i.e., (j, l, k, m, i) = (j, l, k, m, i)
                        temp1 = (gam4d(j,k,m,i)-gam4d(l,k,m,i))/(D(j)-D(l));
                    end %IF ! l
                    dev4(i,j) = dev4(i,j) + L(i,m)*L(m,k)*L(k,l)*L(l,j)*temp1;
                end %DO ! m
            end %DO ! k
        end %DO ! l
    end %DO ! j
end %DO ! i

deri1 = (-(P * gam_tilde_1)/P)*ones(N,1);
deri2 = (2*P* (drv)/P)*ones(N,1);
deri3 = (-6*P* (dev3)/P)*ones(N,1);
deri4 = (24*P* (dev4)/P)*ones(N,1);
deri = [deri1(N/2); deri2(N/2); deri3(N/2); deri4(N/2)];
end