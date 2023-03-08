%Starts the clock to determine the time it takes to run the code%
tic
%clears old variables and the code from the console%
clc;
clear;
format long
saveb1 = [];
savepoles =[];
saveb2 = [];
outhouse =[];
outhome =[];
outsig = [];
%Parameters%
%-----------------------------------------------
winLengthinms=20;
shiftLengthinms=10;
lp_order=20;
mcadams=0.8;
%-----------------------------------------------

%name of the audio file to be processed
filename = 'in.wav';

%this function reads the audio file from the same directory in which the
%matlab script is present. Here, 'sig' stores the audio samples and 'fs' stores
%the sample rate.
[sig,fs] = audioread(filename);

%Incase the audio file has multiple channels, the following code selects
%the first channel only
sig=sig(:,1);

%Improves accuracy by adding 'eps' which is the distance from 1.0 to the
%next largest double precision number.
sig= sig+eps; 

%simulation parameters%
%-----------------------------------------------
%converts the window length from ms to s and multiplies by frame rate to
%get the number of samples in one window
winlen = floor(winLengthinms * 0.001 *fs);
%converts the shift length from ms to s and multiplies by frame rate to
%get the number of samples in one shift window
shift = floor(shiftLengthinms * 0.001 *fs);
%Determines the length of the signal(audio file)
length_sig = numel(sig);
%duraion of the signal
duration_sig = length_sig / fs;
%-----------------------------------------------

%fft processing parameters%
%-----------------------------------------------
NFFT = 2 ^ (ceil(log2(winlen)));
%anaysis and synth window which satisfies the constraint
wPR = hann(winlen);
K = sum(wPR)/shift;
win = sqrt(wPR/K);
Nframes = 1 + floor((length_sig - winlen)/shift);
%carry out the overlap - add FFT processing
sig_rec = zeros(length_sig,1);
%-----------------------------------------------

%Each frame is processed independently, for this a for loop is used.
for m = 1:Nframes-1
    
    
    %indexes for the mth frame are calcualated by the following code:
    %-----------------------------------------------
     concatenated=cat(2,(m*shift+winlen),length_sig);
     adeel = min(concatenated,[],2);
     index = transpose((m*shift):(min(concatenated,[],2))-1); 
     index=index+1;
    %-----------------------------------------------
     %windowed mth frame (other than rectangular window)
     frame = sig(index).*win;
     adeel = frame;
     %linear predictive coding method is used to predict the lpc
     %coeffeicients of the order defined by lp_order. The result is saved
     %into a_lpc.
     [a_lpc,~] = lpc(frame + eps,lp_order);
     saveb1 = [saveb1; a_lpc];
     %The lpc coefficients are converted into poles with the following
     %function.
     [~,poles,~] = tf2zpk(ones(1,numel(a_lpc)),a_lpc);
     savepoles = [savepoles; poles'];
     %imaginary poles are found with the following line of code
     ind_imag = find(imag(poles)~=0 );
     %Only one of the conjugate poles is taken.
     ind_imag_con = ind_imag(1:2:numel(ind_imag)-1);
     %The angle between the positive real axis and the pole is changed by
     %raising the previous angle to the power of the mcadams coefficient.
     %This is done for all the complex poles. The resulting angle defines
     %the new set of poles. 
     new_angles = angle(poles(ind_imag_con)) .^mcadams;
     %A contrainst is set to make the angle between 0 and pi
     new_angles(new_angles(new_angles >= pi)) = pi;
     new_angles(new_angles(new_angles <= 0)) = 0;
     %copy of the original poles to be adjusted with the new angle
     new_poles = poles;
     
     for k = 1:numel(ind_imag_con)-1
         %compute new poles with the same magnitued and new angles
         new_poles(ind_imag_con(k)) = abs(poles(ind_imag_con(k))) * exp(1i*new_angles(k));
         %applied also to the conjugate pole
         new_poles(ind_imag_con(k)+1) = abs(poles(ind_imag_con(k)+1)) * exp(-1i*new_angles(k));
     end
    
     %recover new, modified lpc coefficients
     a_lpc_new = real(poly(new_poles)); 
     %get residual excitation for reconstruction
     res = filter(a_lpc,1,frame);
     %reconstruct frames with new lpc coefficient
     frame_rec = filter(1,a_lpc_new,res);
     saveb2 = [saveb2; frame_rec'];
     frame_rec = frame_rec .* win;
     frame_rec_size = frame_rec;
     outindex = m*shift:((m*shift)+numel(frame_rec))-1;
     outhouse = [outhouse; outindex];
     %Standard overlap add
     sig_rec(outindex)=sig_rec(outindex)+frame_rec;
     outsig = [outsig; sig_rec'];
     outhome = [outhome' sig_rec(outindex)]';
end
maxabs = max(abs(sig_rec));
sig_rec = sig_rec / max(abs(sig_rec));
sig_rec_row = sig_rec';
%The anonymised audio file is written with the changed pole positions. 
audiowrite('out_serial.wav',(sig_rec),fs)
%ends the clock to determine the time it took to run the code%
%fileID = fopen('lpccoeffnotation.txt','w');
%F=[repmat(' %d',1,10),'\n']
%fprintf(fileID,F,saveb1);
%save('lpccoeff.txt','saveb1','-ascii');
writematrix(saveb1,'newlpc.txt');
writematrix(savepoles,'poles.txt');
writematrix(saveb2,'frame_rec.txt');
writematrix(outhome,'outhome.txt');
toc