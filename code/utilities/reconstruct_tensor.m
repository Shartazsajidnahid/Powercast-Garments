function re_myTensor = reconstruct_tensor(M, n_d, r)

re_myTensor = zeros(n_d,size(M.U{2},1),size(M.U{3},1));
for i = 1: n_d
%     if i == n_d - 1
%         pause
%     end
    for j = 1: size(M.U{2},1)
        for k = 1: size(M.U{3},1)
            for r_idx = 1: r
                
                re_myTensor(i,j,k) = re_myTensor(i,j,k) + M.lambda(r_idx) * M.U{1}(i,r_idx) * M.U{2}(j,r_idx) * M.U{3}(k,r_idx);
            end
        end
    end
end