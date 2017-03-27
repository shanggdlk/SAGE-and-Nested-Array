function ret_alpha = compute_alpha(tau, phi, X)

global M F FREQUENCIES L

C = compute_C(repmat(phi, 1,L));
C = squeeze(conj(C(:, 1, :)));

Z_FREQUENCY = repmat(FREQUENCIES, M, 1);
Z_TAU = repmat(tau,M, F);

Z = sum(sum(C .* X .* exp(1j*2*pi.*Z_FREQUENCY.*Z_TAU),1),2);
ret_alpha = 1/(M*F)*Z;

end