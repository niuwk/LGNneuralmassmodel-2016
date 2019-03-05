%% THIS CODE WAS USED IN PART TO GENERATE RESULTS PRESENTED IN THE SPRINGER BOOK CHAPTER, 2016:
% "A Neural Mass Computational Framework to Study Synaptic Mechanisms Underlying Alpha and Theta
% Rhythms", BASABDATTA SEN-BHATTACHARYA and SIMON DURRANT, IN: COMPUTATIONAL NEUROLOGY AND PSYCHIATRY, SPRINGER, 2016
% PLEASE CITE THE WORK IF USING THE CODE IN PART OR FULL.



tic
clear all
% close all

%% MEMBRANE CAPACITANCE IN MICROFARADS
Cm=1;

%% THE TIME PARAMETERS AND VECTOR
delt=0.001; %% 1 millisecond
% endtime=600; %% seconds
endtime=40; %% seconds

timevec=0:delt:endtime;
timelen=numel(timevec);

%% THE STEPS OF SOLUTION
mu=delt;

%% GENERATING NOISY RETINAL INPUT
V_ret_mean=-65;
V_ret_std=2;

% %% TRANSMITTER CONCENTRATION PARAMETERS INITIALISED
Tmax=1; %% in Moles (M): so 2.84mM; from fig. 5 legend 1994 destexhe
Kp = 3.7; %% in mV
Vp=-32; %% in mV
Tval=[Tmax Kp Vp];
%% AMPA PARAMETERS INITIALISED
E_ampa=0;
% scale1=1000;
scale2=50;
alpha_ampa=20*scale2; %% mM^(-1)ms^(-1)------<< BASAL VALUE OS 2: BEING VARIED
beta_ampa=1*scale2;%% ms^(-1)%%
g_ampa_ret2tcr=300; %0.3
g_ampa_ret2in=100; %0.2
g_ampa_tcr2trn=100; %0.2
ampaval=[E_ampa alpha_ampa beta_ampa g_ampa_tcr2trn g_ampa_ret2tcr g_ampa_ret2in];

%% GABA_A PARAMETERS INITIALISED
alpha_gaba_a=20*scale2; %% mM^(-1)ms^(-1)
beta_gaba_a=0.8*scale2;%% ms^(-1)
g_gaba_a_trn2tcr=100;
g_gaba_a_trn2trn=100;
E_gaba_a_trn2tcr=-85; %%
E_gaba_a_trn2trn=-75;
g_gaba_a_in2tcr=g_gaba_a_trn2tcr;
g_gaba_a_in2in=100;%g_gaba_a_trn2trn;
E_gaba_a_in2tcr=E_gaba_a_trn2tcr; %%
E_gaba_a_in2in=E_gaba_a_trn2trn;
gabaAval=[alpha_gaba_a beta_gaba_a g_gaba_a_trn2tcr g_gaba_a_trn2trn E_gaba_a_trn2tcr E_gaba_a_trn2trn g_gaba_a_in2tcr g_gaba_a_in2in E_gaba_a_in2tcr E_gaba_a_in2in];

%% GABA_B PARAMETERS INITIALISED
alpha1_gaba_b=0.2*scale2; %% mM^(-1)ms^(-1)
beta1_gaba_b=0.5*scale2;%% ms^(-1)
alpha2_gaba_b=0.3*scale2; %% mM^(-1)ms^(-1)
beta2_gaba_b=0.1*scale2;%% ms^(-1)
g_gaba_b =60;%% in milliSiemens (mS)

E_gaba_b=-100;
Kd_gaba_b=100;
n=4;

gabaBval=[alpha1_gaba_b  beta1_gaba_b  alpha2_gaba_b  beta2_gaba_b   g_gaba_b   E_gaba_b   Kd_gaba_b   n];
%% LEAK PARAMETERS INITIALISED
E_leak_tcr=-55; % Leak voltage
g_leak_tcr=10; %%

E_leak_trn=-72.5; % Leak voltage
g_leak_trn=10;
% g_leak_trn_arr=[0.0025 0.0075 0.025 0.075 0.25 0.75];

E_leak_in=-72.5; % Leak voltage
g_leak_in=10;

leakval=[E_leak_tcr g_leak_tcr E_leak_trn g_leak_trn E_leak_in g_leak_in];
%% VECTOR INITIALISATION VALUES
% Ginit=1;
Ginit=0.001;

r_ampa_initval=Ginit;
r_gaba_a_initval=Ginit;
r_gaba_b_initval=Ginit;

vinit_tcr=-65; %% in mV
vinit_trn=-85; %% in mV
vinit_in=-75; %% in mV

initval=[Ginit vinit_tcr vinit_trn vinit_in];

% %% CONNECTIVITY PARAMETERS INITIALISED

Cnte=35;
Ctnia=(3/8)*30.9;%% Half of the the inhibitory synaptic connections from the RE pop to the TCR pop
Ctre=7.1;
Cnsi=20;
Cire = 47.4;
Cisi = 23.6;
%% THIS IS THE PARAMETER THAT WE WILL BE DELETING
Ctii=0;%(5/8)*30.9;
Ctnib=0;
Carr=[Cnte Ctnia Ctnib Ctre Cnsi Cire Ctii Cisi];

%% FREQUENCY ANALYSIS
Fs = 1000;%250;
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

%%  MEMBRANE POTENTIAL: TCR, TRN, IN

noftrials=20;
startind=10001; endind=39901;

% EXTRA LOOP ADDED TO TEST WITH VARYING PARAMETERS
%%%%************************************
% VpArr=[-30 -31 -32 -33 -34 -35 -36];
% KpArr=[3 3.2 3.4 3.6 3.8 4];
% 
% for loop = 1:length(VpArr)
%     %% TRANSMITTER CONCENTRATION PARAMETERS INITIALISED
%     Tmax=1; %% in Moles (M): so 2.84mM; from fig. 5 legend 1994 destexhe
%     Kp=KpArr(5);
%     Vp=VpArr(loop);
%     Tval=[Tmax Kp Vp];
    %%%%************************************
    
    Vtcrmat=zeros(noftrials,timelen);
    Vtrnmat=zeros(noftrials,timelen);
    Vinmat=zeros(noftrials,timelen);
    
    for numtrial=1:noftrials
        fprintf('within the %d-th loop \n',numtrial)
        
        %% TRANSMITTER CONCENTRATION VECTOR CORRESPONDING TO RETINAL NOISY INPUT
        R=randn(1,timelen);
        V_ret = ((((R-mean(R)) ./ std(R)) .* V_ret_std) + V_ret_mean);%V_ret_mean.*ones(1,timelen);%
        T_ret=Tmax./(1+exp(-(V_ret-Vp)./Kp));
        
        %% EVENT INPUTS AT DIFFERENT FREQUENCIES
        %       int_hz=167; % interval for 6 Hz
        %       int_hz=125; % interval for 8 Hz
        %       int_hz=100; % interval for 10 Hz
        %       int_hz=62; % interval for 16 Hz
        %       int_hz=40; % interval for 25 Hz
        
        % PIGGY BACKING THE NOISE WITH THE IMPULSE TRAIN
        
        %     Vimpulse = zeros(1, length(V_ret));
        %     Vimpulse(1:int_hz:end)=10;%% NOTE THAT IF THIS IS REDUCED TO 25 (SAY), THERE WILL BE MUCH 'NOISE' IN THE POWER SPECTRUM AND THE HARMONICS WILL NOT BE SO PROMINENT! IS THIS SIMULATING THE CLOSED EYE STIMULATION EXPERIMENT I.E. WITH INPUT OF LESS INTENSITY?
        %     V_eventinp = V_ret + Vimpulse;
        
        
        %     T_ret=Tmax./(1+exp(-(V_eventinp-Vp)./Kp));
        
        
        [Y TrConc]=rk45func_thalmodwithin(Carr, initval, Tval, ampaval, gabaAval, gabaBval, leakval, T_ret, Cm);
        
        Vtcrmat(numtrial,:)=Y(10,1:end-1);
        Vtrnmat(numtrial,:)=Y(11,1:end-1);
        Vinmat(numtrial,:)=Y(12,1:end-1);
        T_tcrArr(numtrial,:)=TrConc(1, 1:end-1);
        T_inArr(numtrial,:)=TrConc(2, 1:end-1);
        T_trnArr(numtrial,:)=TrConc(3, 1:end-1);
        
        
        %% Filter data
        filtData1 = filtfilt(B,A,Vtcrmat(numtrial,startind:endind));
        filtData2 = filtfilt(B,A,Vtrnmat(numtrial,startind:endind));
        filtData3 = filtfilt(B,A,Vinmat(numtrial,startind:endind));
        
        hpopts1 = psdopts(hp,filtData1);
        hpopts2 = psdopts(hp,filtData2);
        hpopts3 = psdopts(hp,filtData3);
        
        set(hpopts1,'Fs',Fs,'NFFT',NFFT,'Normalized',Normalised)
        set(hpopts2,'Fs',Fs,'NFFT',NFFT,'Normalized',Normalised)
        set(hpopts3,'Fs',Fs,'NFFT',NFFT,'Normalized',Normalised)
        
        hpsd1 = psd(hp,filtData1,hpopts1);
        hpsd2 = psd(hp,filtData2,hpopts2);
        hpsd3 = psd(hp,filtData3,hpopts3);
        
        Ptcr(numtrial,:)=hpsd1.Data';
        Ptrn(numtrial,:)=hpsd2.Data';
        Pin(numtrial,:)=hpsd3.Data';
    end
    Vtcravgmat=mean(Vtcrmat,1);
    Vtrnavgmat=mean(Vtrnmat,1);
    Vinavgmat=mean(Vinmat,1);
    Ptcravg=mean(Ptcr,1);
    Ptrnavg=mean(Ptrn,1);
    Pinavg=mean(Pin,1);
    T_mean_mat=zeros(3,length(T_tcrArr));
    T_mean_mat(1,:)=mean(T_tcrArr,1);
    T_mean_mat(2,:)=mean(T_trnArr,1);
    T_mean_mat(3,:)=mean(T_inArr,1);
    toc
    figure, plot(Vtcravgmat(1,startind:endind),T_mean_mat(1,startind:endind),'bo'),
    ylim([-0.001 max(max(T_mean_mat))]); xlim([-90 -35]);
    hold on, plot(Vtrnavgmat(1,startind:endind),T_mean_mat(2,startind:endind),'--g','linewidth',2)
    hold on, plot(Vinavgmat(1,startind:endind),T_mean_mat(3,startind:endind),'r+','markersize',4)
    legend('TCR','TRN','IN')
    xlabel('membrane potential (mV)','Fontsize',14)
    ylabel('Transmitter Concentration (mM)','Fontsize',14)
    axis('xy')
    set(gca,'Fontsize',12),
    
    
    %  figure,plot(V_eventinp)
    
    
    %% Visualising power spectra
    % THE TIME SERIES
    figure, plot(timevec(startind:endind),Vtcravgmat(1,startind:endind))%ylabel('TCR','fontsize',14)
    xlabel('Time (seconds)','fontsize',16)
    axis('xy')
    set(gca,'Fontsize',14),box off
    
    figure, plot(timevec(startind:endind),Vinavgmat(1,startind:endind))%ylabel('TCR','fontsize',14)
    xlabel('Time (seconds)','fontsize',16)
    axis('xy')
    set(gca,'Fontsize',14), box off
    
    figure, plot(timevec(startind:endind),Vtrnavgmat(1,startind:endind)) %ylabel('Membrane Potential (mV)','fontsize',14)
    xlabel('Time (seconds)','fontsize',16)
    axis('xy')
    set(gca,'Fontsize',14),box off
    
    % THE FOURIER TRANSFORM AND POWER SPECTRAL DENSITY
    fr=hpsd1.Frequencies;
    
    figure,plot(fr,Ptcravg),box on, title('TCR')
    xlim([0 50])
    axis('xy')
    set(gca,'Fontsize',12),
    figure,plot(fr,Ptrnavg,'g'),box on, title('TRN')
    xlim([0 50])
    axis('xy')
    set(gca,'Fontsize',12),
    figure,plot(fr,Pinavg,'m'), box on, title('IN')
    %     xlabel('Frequency (Hz)','fontsize',14),
    %     ylabel('Power spectral density','fontsize',14),
    xlim([0 50])
    axis('xy')
    set(gca,'Fontsize',12),
    
%   Ptcravgmat(loop,:)=Ptcravg;
% 
% end
  
    %% SHORT TIME FOURIER TRANSFORM
    % TCR
    Vavgmat=[Vtcravgmat;Vinavgmat;Vtrnavgmat];
    time_window_len=1000;
    slider_span=0.5*time_window_len;
    WIN = transpose(hamming(time_window_len)); %ones(1,time_window_len);
    vismat=[];
    for pass=1:size(Vavgmat,1)
        filtData = filtfilt(B,A,Vavgmat(pass,10000:39000));
        Data=filtData;
        
        x=size(Data);
        no_of_windows=round(x(2)/slider_span);
        
        xind = ceil(x(2)/slider_span);
        yind = (NFFT/2)+1;
        stft_mat = zeros(xind,yind);
        zeropadding_len=time_window_len - slider_span;
        X = [zeros(1,zeropadding_len) Data zeros(1,zeropadding_len)]; % padding the signal with zeros
        % Hamming window eases out the ripples compared to a rectangular window
        iter = 0;
        for i = 1:slider_span:x(2)-1
            iter = iter + 1;
            time_window = X(i:(i + time_window_len - 1));  % making window
            
            signal_window = time_window.* WIN; %   y(t)=h(t)*w(n-t)
            signal_out(iter,:) = abs(fft(signal_window, NFFT)); %fft of y(t) gives us stft
            stft_mat(iter,:) = signal_out(iter,1:yind);
        end
        fr = 0:0.25:Fs/2;
        vismat2=stft_mat';
        
        vismat=[vismat, vismat2];
    end
    vismat = normalise(vismat,100);
    x2=Fc1; y2=Fc2;
    figure, imagesc([],fr((4*x2+1):(4*y2+1)),vismat((4*x2+1):(4*y2+1),2:end));
    xlabel('Time windows','Fontsize',14);
    ylabel('frequency(Hz)','Fontsize',14);
    axis('xy'), colormap('summer')
    set(gca,'Fontsize',12),
    ylim([1 50])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAR CHARTS

delta_low=find(fr==1);
delta_high=find(fr==3.5);
theta_low=find(fr==3.75);
theta_high=find(fr==7.5);
alpha_low=find(fr==7.75);
alpha_high=find(fr==13.5);
beta_low=find(fr==13.75);
beta_high=find(fr==30.5);

    %% TCR

% NORMALISE THE POWER ARRAY FOR TCR BETWEEN 0 AND 1
Ptcravg=normalise(Ptcravg,1);
% NEXT EXTRACT THE AVERAGE POWER IN EACH FREQUENCY BAND
deltarelpow=mean(Ptcravg(1,delta_low:delta_high),2);
thetarelpow=mean(Ptcravg(1,theta_low:theta_high),2);
alpharelpow=mean(Ptcravg(1,alpha_low:alpha_high),2);
betarelpow=mean(Ptcravg(1,beta_low:beta_high),2);
% NOW FORM THE ARRAY FOR MAKING THE BAR CHART CONSISTING OF THE FOUR
% FREQUENCY COMPONENTS
relpowarr_tcr(:,1)=[deltarelpow; thetarelpow; alpharelpow; betarelpow];
% NOW DRAW THE FIGURE
figure,
   bar1=bar(transpose(relpowarr_tcr),'hist');colormap('jet')
   xlim([0.45 4.55]), ylabel('Normalised Power','Fontsize',14)
  
% 
% % 


