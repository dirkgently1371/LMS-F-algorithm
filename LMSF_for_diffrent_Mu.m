
clear
clc
tic
%--------------------------------------------------------------------------
% initial values

% Number of iterations and samples
n = 20000;
% L is length of Filter
L = 16;
% SNR_db : [5 10 15 25]
SNR_dB = 5;
% thereshold parameter (for high SNR we should set epsilon in range of (0:0.3))
% thereshold parameter (for high SNR we should set epsilon in range of (0.3:1))
epsilon = 0.9;
% Number of Experiments
T = 1000;
% Mu is step size parameter for LMSF algorithm
Mu = [1e-2 8e-3 6e-3 4e-3 2e-3 1e-3 8e-4 6e-4 4e-4 2e-4 1e-4];
% System_Vector is coefficient of channel
System_Vector = nan(1,L);
%--------------------------------------------------------------------------
% we assume a constant system

% nn and b are parameters of system
nn = 1:16;
b = 3.5;
for i = 1 : L
    System_Vector(1,i) = .5 * (1 + cos(2 * pi * (i - 2) / b));
end
%--------------------------------------------------------------------------
Wtemp_LMSF = 0;
for z=1:length(Mu)
    
    % t is number of Experiments
    for t = 1:T
        fprintf('%d from %d Experiments for Mu = %f\n',t,T,Mu(z))
        % Input_Signal is input of unknown system : X(n)
        % Input_Signal modulated by BPSK
        Input_Signal = (2 * mod(reshape(randperm(n*1),n,1),2) - 1)';
        System_Output = filter(System_Vector,1,Input_Signal);
        % NoisyOutput = System_Output + Noise
        NoisyOutput = awgn(System_Output, SNR_dB);
        % initialize of LMF parameters
        W_vector_LMSF = zeros(1,L);
        Input_Vector_LMSF = zeros(1,L);
        
        for j = 1 : n
            
            % LMSF Algorithm
            Input_Vector_LMSF(1,2:end) = Input_Vector_LMSF(1,1:end-1);
            Input_Vector_LMSF(1,1) = Input_Signal(j);
            out_LMSF = (W_vector_LMSF) * Input_Vector_LMSF';
            Error_LMSF = NoisyOutput(j) - out_LMSF;
            W_vector_LMSF = W_vector_LMSF +  Mu(z) * (Error_LMSF .^ 3) *...
                Input_Vector_LMSF / ((Error_LMSF .^ 2) + epsilon);
            temp_LMSF(j,t) = norm((W_vector_LMSF-System_Vector),2) ^ 2; %#ok
        end
        
        Wtemp_LMSF = Wtemp_LMSF + W_vector_LMSF;
        
    end
    
    % MSD calculation
    for i=1:n
        
        msd_lmsf(i) = mean(temp_LMSF(i,:) / norm(System_Vector) ^ 2);   %#ok
    end
    
    W_vector_LMSF = Wtemp_LMSF ./ T;
    
    semilogy(msd_lmsf,'LineWidth',1.3);
    hold on
end

%--------------------------------------------------------------------------
% plot part

lgh=legend(strcat('LMS: \mu=5e-2'),strcat('LMF: \mu=2e-3'),strcat('LMS/F \mu=1e-2')...
    ,strcat('LMS/F \mu=8e-3'),strcat('LMS/F \mu=6e-3'),strcat('LMS/F \mu=4e-3')...
    ,strcat('LMS/F \mu=2e-3'),strcat('LMS/F \mu=1e-3'),strcat('LMS/F \mu=8e-4')...
    ,strcat('LMS/F \mu=6e-4'),strcat('LMS/F \mu=4e-4'),strcat('LMS/F \mu=2e-4')...
    ,strcat('LMS/F \mu=1e-4'));
title(['Learning Curve for LMSF,LMF,LMS',10,'(SNR = ',num2str(SNR_dB),...
    ' dB',' , ','Filter-Length = ',num2str(L),...
    ' , ','\epsilon = ',num2str(epsilon),'  ','Number of Experiment = ',num2str(n),' )'])
toc