
clear
clc
tic
%--------------------------------------------------------------------------
% initial values

% Number of iterations and samples
n = 20000;
% Number of Experiments
T = 1000;
% SNR_db : [5 10 15 25]
SNR_dB = 5;
% L is length of Filter
L = 16;
% Mu is step size parameter for LMS algorithms
Mu1 = 0.05;
% Mu is step size parameter for LMF algorithms
Mu2 = 0.002;
% System_Vector is coefficient of channel
System_Vector = nan(1,L);
%--------------------------------------------------------------------------
% we assume a constant system

% nn and b are parameters of system
nn = 1:16;
b = 3.5;
for i = 1:L
    System_Vector(1,i) = .5 * (1 + cos(2 * pi * (i - 2) / b));
end
%--------------------------------------------------------------------------
Wtemp_LMS = 0;
Wtemp_LMF = 0;
% t is number of Experiments
for t = 1:T
    fprintf('%d from %d Experiments\n',t,T)
    % Input_Signal is input of unknown system : X(n)
    % Input_Signal modulated by BPSK
    Input_Signal = (2 * mod(reshape(randperm(n*1),n,1),2) - 1)';
    System_Output = filter(System_Vector,1,Input_Signal);
    % NoisyOutput = System_Output + Noise
    NoisyOutput = awgn(System_Output, SNR_dB);
    % initialize of LMS parameters
    W_vector_LMS = zeros(1,L);
    Input_Vector_LMS = zeros(1,L);
    % initialize of LMF parameters
    W_vector_LMF = zeros(1,L);
    Input_Vector_LMF = zeros(1,L);
    
    for j = 1:n
        
        % LMS Algorithm
        
        Input_Vector_LMS(1,2:end) = Input_Vector_LMS(1,1:end-1);
        Input_Vector_LMS(1,1) = Input_Signal(j);
        out_LMS = (W_vector_LMS) * Input_Vector_LMS';
        Error_LMS = NoisyOutput(j) - out_LMS;
        W_vector_LMS = W_vector_LMS +  Mu1 * Error_LMS * Input_Vector_LMS;
        
        % LMF Algorithm
        Input_Vector_LMF(1,2:end) = Input_Vector_LMF(1,1:end-1);
        Input_Vector_LMF(1,1) = Input_Signal(j);
        out_LMF = (W_vector_LMF)*Input_Vector_LMF';
        Error_LMF = NoisyOutput(j) - out_LMF;
        W_vector_LMF = W_vector_LMF +  Mu2 * (Error_LMF .^ 3) * Input_Vector_LMF;
        temp_LMS(j,t) = norm((W_vector_LMS-System_Vector),2) ^ 2;  %#ok
        temp_LMF(j,t) = norm((W_vector_LMF-System_Vector),2) ^ 2;  %#ok
        
    end
    
    Wtemp_LMS = Wtemp_LMS + W_vector_LMS;
    Wtemp_LMF = Wtemp_LMF + W_vector_LMF;
    
end

% MSD calculation
for i=1:n
    
    msd_lms(i) = mean(temp_LMS(i,:) / norm(System_Vector) ^ 2);    %#ok
    msd_lmf(i) = mean(temp_LMF(i,:) / norm(System_Vector) ^ 2);    %#ok
    
end

W_vector_LMS = Wtemp_LMS ./ T;
W_vector_LMF = Wtemp_LMF ./ T;
%--------------------------------------------------------------------------
% plot part

semilogy(msd_lms,'LineWidth',1.3);
hold on
semilogy(msd_lmf,'LineWidth',1.3);
lgh=legend(strcat('LMS: \mu=5e-2'),strcat('LMF: \mu=2e-3'),strcat('LMS/F \mu=1e-2')...
    ,strcat('LMS/F \mu=8e-3'),strcat('LMS/F \mu=6e-3'),strcat('LMS/F \mu=4e-3')...
    ,strcat('LMS/F \mu=2e-3'),strcat('LMS/F \mu=1e-3'),strcat('LMS/F \mu=8e-4')...
    ,strcat('LMS/F \mu=6e-4'),strcat('LMS/F \mu=4e-4'),strcat('LMS/F \mu=2e-4')...
    ,strcat('LMS/F \mu=1e-4'));
title(['Learning Curve for LMF,LMS',10,'(SNR = ',num2str(SNR_dB),...
    ' dB',' , ','Filter-Length = ',num2str(L),' , ',...
    'Number of Experiment = ',num2str(T),' )'])
grid on
xlabel(' n (iteration)')
ylabel(' MSD ')
toc