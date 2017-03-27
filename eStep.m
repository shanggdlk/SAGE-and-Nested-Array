function X = eStep(parameter, csi, K)

global M L FREQUENCIES F

S_FREQUENCY = repmat(reshape(FREQUENCIES, [1,1,F]), M, L);

    
%% compute the X_k;
alpha = repmat(parameter.alpha, M, 1, F);

C = compute_C(parameter.phi);

S_TAU = repmat(parameter.tau, M, 1, F);

S_matrix = alpha.* C .* exp(-1j*2*pi.*S_FREQUENCY.*S_TAU);

X = csi - squeeze(sum(S_matrix,2)) + squeeze(S_matrix(:,K,:));

end

