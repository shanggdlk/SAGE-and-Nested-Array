function ret_phi = opt_phi(tau, X)

global DOMAIN_PHI M LAMBDAS F FREQUENCIES BASELAYOUT ANTNUM 

C_M = repmat(transpose(BASELAYOUT),1,DOMAIN_PHI.length, F);
[~,S_matrix,~] = size(C_M);
phi_space = repmat((DOMAIN_PHI.start:DOMAIN_PHI.step:DOMAIN_PHI.end), ANTNUM,1, F);
C_L = cos(phi_space);
C_LAMBDA = repmat(reshape(LAMBDAS, [1, 1, F]), ANTNUM, DOMAIN_PHI.length, 1);

C = exp(1j*2*pi./C_LAMBDA.*C_M.*C_L); 

C_temp = zeros(M, S_matrix, F);
for K = 1:F
    h = C(:,:,K);
    C_temp(:,:,K) = kr(conj(h),h);
end


C = conj(C_temp);

X = repmat(reshape(X,[M,1,F]),1,DOMAIN_PHI.length, 1);

Z_FREQUENCY = repmat(reshape(FREQUENCIES, [1,1,F]), M, DOMAIN_PHI.length, 1);
Z_TAU = repmat(tau,M,DOMAIN_PHI.length, F);
Z_abs = abs(squeeze(sum(sum(C .* X .* ...
    exp(1j*2*pi.*Z_FREQUENCY.*Z_TAU),1), 3)));
[~, I] = max(Z_abs);

ret_phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;

end