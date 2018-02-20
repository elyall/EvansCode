% SSN

N = 100; % # of steps
numE = 2; % # E units
numI = 2; % # I units
num = numE + numI; % total # of units
labels = [strcat('E',num2str((1:numE)'));strcat('I',num2str((1:numI)'))];

% Build weight matrix
W_EE = 0.044; % excitatory->excitatory
W_EI = 0.023; % excitatory->inhibitory
W_IE = 0.042; % inhibitory->excitatory
W_II = 0.018; % inhibitory->inhibitory
W_EE = W_EE*ones(numE,numE);
W_EI = W_EI*ones(numI,numE); 
W_IE = W_IE*ones(numE,numI);
W_II = W_II*ones(numI,numI);
W = [W_EE,W_IE;W_EI,W_II]; % build connection matrix
% W = W+rand(num); % add randomness to weights
% W = 5;
% W = W*ones(num); % build connection matrix
% W(logical(eye(size(W,1)))) = 0; % set self-connections to zero
W(:,end/2+1:end) = W(:,end/2+1:end)*-1; % make inhibitory connections negative

% firing rates
r_start = 0; % starting firing rate
r = r_start*ones(num,1); % firing rates
n = 2; % power law
k = .04; % constant

% feed-forward strength
hE = 2; % E input strength
hI = 1; % I input strength
h = [hE*ones(numE,1);hI*ones(numI,1)]; % input strength

% time constants
T_E = .02; % time-constant for E units
T_I = .01; % time-constant for I units
T = diag([T_E*ones(1,numE),T_I*ones(1,numI)]); % time constants

% run sim
out = nan(num,N);
out(:,1) = r;
for t = 1:N-1
    x = W*r+h;
    x(x<0) = 0;
    dr = inv(T)*k*(x.^n);
    r = r + dr;
    out(:,t+1) = r;
end
figure;
plot(out');
legend(labels);
