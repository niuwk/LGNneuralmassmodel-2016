%% LOG OCTOBER 2016: 4TH REVISION OF PAPER - AND REVISITING THE RESULTS
% CORRESPONDING TO VARYING CONNECTIVITY PARAMETERS - TABLE 3 IN PAPER. MADE
% A FEW CHANGES TO VISUALISATION AND GENERATION OF DATA: (1) SIMULATION
% TIME IS 40 SECONDS WITH CUT OFF AT 10 AND 30 SECONDS.

% THE CODE IS COPYRIGHT OF BASABDATTA SEN BHATTACHARYA; PLEASE ACKNOWLEDGE THE
% AUTHOR IF USING THE CODE IN PART OR FULL.


tic
clear all
% close all

%% MEMBRANE CAPACITANCE IN MICROFARADS
Cm=1;

%% THE TIME PARAMETERS AND VECTOR
delt=0.001; %% 1 millisecond

endtime=40; %% seconds

timevec=0:delt:endtime;
timelen=numel(timevec);

%% THE STEPS OF SOLUTION
mu=delt;

%% GENERATING NOISY RETINAL INPUT
V_ret_mean=-65; % mV
V_ret_std=2; % mV

%% TRANSMITTER CONCENTRATION PARAMETERS INITIALISED
Tmax=1; %% milliMoles (mM)
Kp = 3.8; %%  mV
Vp=-32; %% mV
Tval=[Tmax Kp Vp];
%% AMPA PARAMETERS INITIALISED
E_ampa=0; %% AMPA reverse potential, mV

alpha_ampa=1000; %% mM^(-1)ms^(-1)
beta_ampa=50;%% ms^(-1)%%
g_ampa_ret2tcr=300; % micro-Siemens/cm^2
g_ampa_ret2in=100; % micro-Siemens/cm^2
g_ampa_tcr2trn=100; % micro-Siemens/cm^2
ampaval=[E_ampa alpha_ampa beta_ampa g_ampa_tcr2trn g_ampa_ret2tcr g_ampa_ret2in];

%% GABA_A PARAMETERS INITIALISED
alpha_gaba_a=1000; %% mM^(-1)ms^(-1)
beta_gaba_a=40;%% ms^(-1)
g_gaba_a_trn2tcr=100; % micro-Siemens/cm^2
g_gaba_a_trn2trn=100;% micro-Siemens/cm^2
E_gaba_a_trn2tcr=-85; %  mV
E_gaba_a_trn2trn=-75; % mV
g_gaba_a_in2tcr=g_gaba_a_trn2tcr; % micro-Siemens/cm^2
g_gaba_a_in2in=100; % micro-Siemens/cm^2
E_gaba_a_in2tcr=E_gaba_a_trn2tcr; %mV
E_gaba_a_in2in=E_gaba_a_trn2trn; %mV
gabaAval=[alpha_gaba_a beta_gaba_a g_gaba_a_trn2tcr g_gaba_a_trn2trn E_gaba_a_trn2tcr E_gaba_a_trn2trn g_gaba_a_in2tcr g_gaba_a_in2in E_gaba_a_in2tcr E_gaba_a_in2in];

%% GABA_B PARAMETERS INITIALISED
alpha1_gaba_b=10; %% mM^(-1)ms^(-1)
beta1_gaba_b=25;%% ms^(-1)
alpha2_gaba_b=15; %% mM^(-1)ms^(-1)
beta2_gaba_b=5;%% ms^(-1)
g_gaba_b =60;% micro-Siemens/cm^2

E_gaba_b=-100; % mV
Kd_gaba_b=100;
n=4;

gabaBval=[alpha1_gaba_b  beta1_gaba_b  alpha2_gaba_b  beta2_gaba_b   g_gaba_b   E_gaba_b   Kd_gaba_b   n];
%% LEAK PARAMETERS INITIALISED
E_leak_tcr=-55; % Leak voltage, mV
g_leak_tcr=10; % micro-Siemens/cm^2

E_leak_trn=-72.5; % Leak voltage, mV
g_leak_trn=10; % micro-Siemens/cm^2
% g_leak_trn_arr=[0.0025 0.0075 0.025 0.075 0.25 0.75];

E_leak_in=-72.5; % Leak voltage, mV
g_leak_in=10; % micro-Siemens/cm^2

leakval=[E_leak_tcr g_leak_tcr E_leak_trn g_leak_trn E_leak_in g_leak_in];
%% VECTOR INITIALISATION VALUES

Ginit=0.001; % Intial value for solving differential equations

r_ampa_initval=Ginit;
r_gaba_a_initval=Ginit;
r_gaba_b_initval=Ginit;

vinit_tcr=-65; %% mV
vinit_trn=-85; %% mV
vinit_in=-75; %% mV

initval=[Ginit vinit_tcr vinit_trn vinit_in];

% %% CONNECTIVITY PARAMETERS INITIALISED

Cnte=35;
Ctnia=(3/8)*30.9;%% Half of the the inhibitory synaptic connections from the RE pop to the TCR pop
Ctnib=(1/8)*30.9;
Ctre=7.1;
Cnsi=20;
Cire = 47.4;
Cisi =23.6;
%% THIS IS THE PARAMETER THAT WE WILL BE SETTING TO ZERO TO TEST REMOVAL OF IN FROM THE LGN NETWORK
Ctii = 15.45;

Carr=[Cnte Ctnia Ctnib Ctre Cnsi Cire Ctii Cisi];

%% EVOLUTION OF THE MEMBRANE POTENTIALS

noftrials=20;
startind=9001; endind=39001;

Vtcrmat=zeros(noftrials,timelen);
Vtrnmat=zeros(noftrials,timelen);
Vinmat=zeros(noftrials,timelen);

%% FREQUENCY ANALYSIS
Fs = 1000;%Sampling frequency
NFFT=4*Fs;
WindowType = 'hamming';
SegmentLength=(1/2)*Fs;
OverlapPercent=50;
Normalised=0;
hp = spectrum.welch(WindowType,SegmentLength,OverlapPercent);


%% Construct an FDESIGN object and call its BUTTER method.
Fc1 = 1; % First Cutoff Frequency
Fc2 = 100; % Second Cutoff Frequency
N = 10; % Order
h = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');
[B,A]=sos2tf(Hd.sosMatrix,Hd.Scalevalues);


for numtrial=1:noftrials
    fprintf('within the %d-th loop \n',numtrial)
    
    %% TRANSMITTER CONCENTRATION VECTOR CORRESPONDING TO RETINAL NOISY INPUT
    R=randn(1,timelen);
    V_ret = ((((R-mean(R)) ./ std(R)) .* V_ret_std) + V_ret_mean);
    
     
    T_ret=Tmax./(1+exp(-(V_ret-Vp)./Kp));
    
    
    Y=rk45func_thalmodwithin(Carr, initval, Tval, ampaval, gabaAval, gabaBval, leakval, T_ret, Cm);
    
    Vtcrmat(numtrial,:)=Y(10,1:end-1);%% TCR
    Vtrnmat(numtrial,:)=Y(11,1:end-1);%% TRN
    Vinmat(numtrial,:)=Y(12,1:end-1);%% IN
    
    
    
    %% Filter data
    filtData1 = filtfilt(B,A,Vtcrmat(numtrial,startind:endind));%% TCR
    filtData2 = filtfilt(B,A,Vtrnmat(numtrial,startind:endind));%% TRN
    filtData3 = filtfilt(B,A,Vinmat(numtrial,startind:endind));%% IN
    
    %% Use filtered data to calculate power spectrum
    hpopts1 = psdopts(hp,filtData1);%% TCR
    hpopts2 = psdopts(hp,filtData2);%% TRN
    hpopts3 = psdopts(hp,filtData3);%% IN
    
    set(hpopts1,'Fs',Fs,'NFFT',NFFT,'Normalized',Normalised)%% TCR
    set(hpopts2,'Fs',Fs,'NFFT',NFFT,'Normalized',Normalised)%% TRN
    set(hpopts3,'Fs',Fs,'NFFT',NFFT,'Normalized',Normalised)%% IN
    
    hpsd1 = psd(hp,filtData1,hpopts1);%% TCR
    hpsd2 = psd(hp,filtData2,hpopts2);%% TRN
    hpsd3 = psd(hp,filtData3,hpopts3);%% IN
    
    Ptcr(numtrial,:)=hpsd1.Data';%% TCR
    Ptrn(numtrial,:)=hpsd2.Data';%% TRN
    Pin(numtrial,:)=hpsd3.Data';%% IN
end
%% Calculate the mean membrane potential across trials
Vtcravgmat=mean(Vtcrmat,1);%% TCR
Vtrnavgmat=mean(Vtrnmat,1);%% TRN
Vinavgmat=mean(Vinmat,1);%% IN
%% Calculate the mean power spectrum across trials
Ptcravg=mean(Ptcr,1);%% TCR
Ptrnavg=mean(Ptrn,1);%% TRN
Pinavg=mean(Pin,1);%% IN
toc


%% Visualising the data
% THE TIME SERIES
figure,subplot(3,1,1), plot(timevec(startind:endind),Vtcravgmat(1,startind:endind),'r'), title('TCR')
xlabel('Time (seconds)','fontsize',14)

set(gca,'Fontsize',12),box off

hold on, subplot(3,1,2), plot(timevec(startind:endind),Vinavgmat(1,startind:endind),'m'), title('IN')
xlabel('Time (seconds)','fontsize',14)

set(gca,'Fontsize',12), box off

hold on, subplot(3,1,3), plot(timevec(startind:endind),Vtrnavgmat(1,startind:endind),'b'), title('TRN')
xlabel('Time (seconds)','fontsize',14)

set(gca,'Fontsize',12),box off

% THE FOURIER TRANSFORM AND POWER SPECTRAL DENSITY
fr=hpsd1.Frequencies;

figure, subplot (3,1,1), plot(fr,Ptcravg,'r'),box on, title('TCR')
xlim([0 50])

set(gca,'Fontsize',12),

hold on, subplot(3,1,2), plot(fr,Pinavg,'m'), box on, title('IN')
xlim([0 50])

set(gca,'Fontsize',12),

hold on, subplot(3,1,3), plot(fr,Ptrnavg,'b'),box on, title('TRN')
xlim([0 50])
axis('xy')
set(gca,'Fontsize',12)