function [parameters,ret_parameter] = channelEstimation(csi_sample)
%% csi_sample: M*1 complex
% theta is a L*3 vector: [alpha, phi, tau] 
% [alpha, -, -]_l: the amplitude of the path l;
% [-, phi, -]_l: the azimuth angle of the path l;
% [-, -, tau]_l: the delay of the path l;

globals_init();

%  no input csi -> simulation
if nargin <= 0
    [csi_sample, ~] = generate_simulation();
end



% disp(csi_sample);

global ITERATION L DOMAIN_TAU M F SIMULATION_PHI SIMULATION_TAU SPEED_OF_LIGHT...

%disp(csi_sample);
% change the CSI;
new_csi = zeros(M,F);
for I =  1:F
    h = csi_sample(:,I);
    s = reshape(transpose(h*h'),[],1);
    new_csi(:,I) = s;
end

csi_sample = new_csi;
% disp(csi_sample);
% G = unique(abs(roundn(csi_sample,-3)), 'rows');
% disp(G);

% parameters: ITERATION+1 cell, inside each cell is a struct 
% (alpha, phi, tau)
parameters = cell(ITERATION+1, 1);

for i = 1:ITERATION+1
    parameters{i} = struct('alpha', zeros(1, L), ...
    'phi', zeros(1, L), 'tau', zeros(1, L));
end

% parameters{1}.tau = zeros(1, L)+...
%     DOMAIN_TAU.start+DOMAIN_TAU.step*round(DOMAIN_TAU.length/2);
parameters{1} = init(csi_sample, parameters{1});
% groundtruth as input
% parameters{1}.tau = SIMULATION_TAU;
% parameters{1}.phi = SIMULATION_PHI;
% parameters{1}.alpha = (1+1j)./(SPEED_OF_LIGHT*SIMULATION_TAU);

%% Iterating
for I = 1:ITERATION
    for K = 1:L
    % Expectation step (E-step)
    X = eStep(parameters{I}, csi_sample, K);
    %disp(X);
    % Maximization step (M-step)
    parameters{I+1}.tau(K) = opt_tau(parameters{I}.phi(K), X);
    parameters{I+1}.phi(K) = opt_phi(parameters{I+1}.tau(K), X);
    parameters{I+1}.alpha(K) = compute_alpha(parameters{I+1}.tau(K),...
        parameters{I+1}.phi(K), X);
    
    parameters{I}.tau(K) = parameters{I+1}.tau(K);
    parameters{I}.phi(K) = parameters{I+1}.phi(K);
    parameters{I}.alpha(K) = parameters{I+1}.alpha(K);
    end

    %% compute expectation
%     alpha = parameters{I+1}.alpha(K);
%     alpha = repmat(alpha,M,1);
% 
%     phi = parameters{I+1}.phi(K);
%     phi = repmat(phi,M,1);
% 
%     tau = parameters{I+1}.tau(K);
%     tau = repmat(tau,M,1);
% 
%     C_M = repmat(transpose(0:M-1), 1, L);
%     C = exp(1j*2*pi/LAMBDA*D*C_M.*cos(phi));
% 
%     S = sum(alpha.*C.*exp(-1j*2*pi*tau*FREQUENCY),2);
    %disp(csi_sample-S);
end

ret_parameter = parameters{ITERATION};
disp(ret_parameter.tau);
disp(ret_parameter.phi);
%disp(ret_parameter.alpha);
disp(SIMULATION_TAU);
disp(SIMULATION_PHI);

% plot(h(1,:));
% hold on;
% plot(h(2,:));
% plot(h(3,:));
% plot(h(4,:));
% hold off;


% H = [];
% for G = 1:ITERATION/L  
%     H = [H;parameters{G*L}.tau];
% end
% plot(H);

end


%%%%%%%%%%%%%%%%%%%%%code stop here

function [csi_sample,csi_hidden] = generate_simulation()
    global SIMULATION_TAU SIMULATION_PHI BASELAYOUT L SPEED_OF_LIGHT ANTNUM FREQUENCIES F LAMBDAS
    
    C_M = repmat(transpose(BASELAYOUT), 1, L, F);
    C_L = repmat(cos(SIMULATION_PHI), ANTNUM, 1, F);
    C_LAMBDA = repmat(reshape(LAMBDAS, [1,1,F]), ANTNUM, L, 1);
    %disp(size(C_L));
    C = exp(1j*2*pi./C_LAMBDA.*C_M.*C_L);
    
    
    %disp(size((1+1j)./(SPEED_OF_LIGHT*SIMULATION_TAU)));
    ALPHA = repmat((1+1j)./(SPEED_OF_LIGHT*SIMULATION_TAU), ANTNUM, 1, F);
    CSI_TAU = repmat(SIMULATION_TAU, ANTNUM, 1, F);
    CSI_FREQUENCY = repmat(reshape(FREQUENCIES, [1,1,F]), ANTNUM, L);
    %disp(size(ALPHA));
    csi_sample = ALPHA .* C .* exp(-1j*2*pi.*CSI_TAU.*CSI_FREQUENCY);
    %disp(C);
    csi_hidden = csi_sample;
    csi_sample = squeeze(sum(csi_sample, 2));
    %disp(csi_sample);
end




function globals_init
%% physical
% LAMBDA: wavelength of the signal;
% FREQUENCY: frequency of the signal;
% M: antenna array size;
% L: # propagation paths;
% D: spacing between adjacent rx antenna
% N: # sample
% F: # measured subcarrier
% DELTA_FREQUENCY: difference between adjacent subcarrier
    global CENTRAL_FREQUENCY SPEED_OF_LIGHT M L D ITERATION DOMAIN_TAU ...
        DOMAIN_PHI F DELTA_FREQUENCY FREQUENCIES LAMBDAS ANTNUM BASELAYOUT
    CENTRAL_FREQUENCY = 5.2e9;  %unit hz
    DELTA_FREQUENCY = 20e6/64; % 20Mhz, 64 slots
    F = 56;    
    SPEED_OF_LIGHT = 3e8;  %unit m/s
    
    FREQUENCIES = ((1:F) - (F+1)/2) * DELTA_FREQUENCY+CENTRAL_FREQUENCY;
    LAMBDAS = SPEED_OF_LIGHT./FREQUENCIES;

    ANTNUM = 4;
    M = ANTNUM^2;
    L = 4;
    D = mean(LAMBDAS)/2;
    ITERATION = 50;
    DOMAIN_TAU = struct('start', 10e-9, 'end', 30e-9, 'step', 1e-9); % unit: s
    DOMAIN_TAU.length = round((DOMAIN_TAU.end - DOMAIN_TAU.start) ...
        / DOMAIN_TAU.step + 1);
    
    DOMAIN_PHI = struct('start', 0, 'end', pi, 'step', pi/100); % unit: radius
    DOMAIN_PHI.length = round((DOMAIN_PHI.end - DOMAIN_PHI.start) ...
        / DOMAIN_PHI.step + 1);
    BASELAYOUT = D*[0,1,2,5];
    %BASELAYOUT = D*[0, 1, 2, 3, 7, 11];
    %BASELAYOUT = D*[0,1,2,3,4,9,14,19];
    %BASELAYOUT = D*[0,3,6,9,5,10,15,20];
    %BASELAYOUT = D*[0,1,2,3,4,5,6,7];
    %% simulation
    global SIMULATION_TAU SIMULATION_PHI
    %SIMULATION_TAU = [12 13 15 17 19 20 22 24]*1e-9;
    %SIMULATION_PHI = [0.1 0.3 0.4 0.45 0.5 0.6 0.76 0.8]*pi;
    SIMULATION_TAU = [12 13 15 17]*1e-9;
    SIMULATION_PHI = [0.1 0.3 0.4 0.6]*pi;
    
    %
   
end