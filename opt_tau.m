function ret_tau = opt_tau(phi, X)

global DOMAIN_TAU M FREQUENCIES L F

% only need a slice of C
C = compute_C(repmat(phi, 1, L));
C = conj(C(:, 1, :));

C = repmat(C,1,DOMAIN_TAU.length, 1);

X = repmat(reshape(X, [M, 1, F]),1,DOMAIN_TAU.length, 1);

Z_FREQUENCY = repmat(reshape(FREQUENCIES, [1,1,F]), M, DOMAIN_TAU.length);
Z_TAU = repmat(DOMAIN_TAU.start:DOMAIN_TAU.step:DOMAIN_TAU.end,...
    M,1, F);

Z_abs = squeeze(abs(sum(sum(C.*X.*exp(1j*2*pi.*Z_FREQUENCY.*Z_TAU),1),3)));
[~, I] = max(Z_abs);

ret_tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

end