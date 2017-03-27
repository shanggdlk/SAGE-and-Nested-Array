function ret_tau = init_tau(X)

global DOMAIN_TAU M FREQUENCIES F 

X = repmat(reshape(X, [M, 1, F]),1,DOMAIN_TAU.length, 1);

FREQUENCY_matrix = repmat(reshape(FREQUENCIES, [1,1,F]), M, DOMAIN_TAU.length, 1);
TAU_matrix = repmat(DOMAIN_TAU.start:DOMAIN_TAU.step:DOMAIN_TAU.end,...
    M,1, F);

%Tau_val = squeeze(sum(sum(abs(X.*exp(1j*2*pi.*FREQUENCY_matrix.*TAU_matrix)),1),3));
Tau_val = squeeze(sum(abs(sum(X.*exp(1j*2*pi.*FREQUENCY_matrix.*TAU_matrix),3)),1));
[~, I] = max(Tau_val);

ret_tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

end