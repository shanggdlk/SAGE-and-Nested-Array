function ret = compute_C(phi)
%M * L * F
global L LAMBDAS F BASELAYOUT ANTNUM M

C_M = repmat(transpose(BASELAYOUT), 1, L, F);
C_L = repmat(cos(phi), ANTNUM, 1, F);
C_LAMBDA = repmat(reshape(LAMBDAS, [1,1,F]), ANTNUM, L, 1);
C_temp = exp(1j*2*pi./C_LAMBDA.*C_M.*C_L);

C = zeros(M, L, F);

for K = 1:F
    h = C_temp(:,:,K);
    C(:,:,K) = kr(conj(h),h);
end

ret = C;

end