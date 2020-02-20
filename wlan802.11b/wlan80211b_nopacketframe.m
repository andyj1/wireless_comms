% ECE408 - Wireless Communications
% Jongoh (Andy) Jeong
% 802.11b WLAN Standard - Simulation Project
% Date: February 19, 2020
clear all; close all; clc;

%% IEEE 802.11b
reset(RandStream.getGlobalStream);
rng default; % for reproducibility

% Simulation Parameters
snrVector = -8:2:10;
nIter = 1; 
% min: 4, max: 8192 as per 802.11 standards
PacketSizeBits = 8192; % 802.11 packet size
% samples per chip; samples will be up/down-sampled to this reference
spc = 8;       

% call objects for modulation / demodulation with appropriate parameters
modFcns = {@(x, rate) ModSchemes.BarkerModulator(x,rate), ...
           @(x, rate) ModSchemes.BarkerModulator(x,rate), ...
           @(x, rate) ModSchemes.CCKModulator(x,rate), ...
           @(x, rate) ModSchemes.CCKModulator(x,rate)};
demodFcns = {@(x, rate) ModSchemes.BarkerDemodulator(x,rate), ...
             @(x, rate) ModSchemes.BarkerDemodulator(x,rate), ...
             @(x, rate) ModSchemes.CCKDemodulator(x,rate), ...
             @(x, rate) ModSchemes.CCKDemodulator(x,rate)};

% supported data rates as per 802.11b standard
dataRates = [1, 2, 5.5, 11];
% bits per symbol for all data rates
BPSes = [1, 2, 4, 8];
% chip spreading rates for all data rates
chipSpreadLengths = [11, 11, 8, 8];

BERVector = zeros(length(dataRates),length(snrVector));
fprintf('Simulation starting...\n'); tic;
for rate=1:length(dataRates) %corrersponding to 4 different rate options (1, 2, 5.5, 11 Mbps)
    fprintf('Data rate: %.f\n', dataRates(rate));
    for i=1:length(snrVector)
        fprintf('SNR: %f', snrVector(i));
        totalbits = 0; 
        nerror = 0;
        for iter = 1:nIter
            % adjust SNR
            sampRate = chipSpreadLengths(rate) * spc;
            snrdB = snrVector(i) + 10*log10(BPSes(rate)) - 10*log10(sampRate);
            % generate packet
            msgBin = randi([0 1],PacketSizeBits,1);
            % modulate
            msgMod = modFcns{rate}(msgBin, dataRates(rate));
            % upsample, pass through a pulse shaping filter
            [h, upsampledChips, chipFilterDelay] = Filter.PulseShapeFilter(msgMod, spc); 
            samples = filter(h,1,upsampledChips);
            % pass through an AWGN channel with adjusted SNR
            txNoisy = awgn(samples, snrdB, 'measured');
            % filer back (and downsample), perfectly knowing impulse response of the filter
            filtTxSig = filter(h,1,txNoisy); % column vector
            [RxChips,bitDelay] = Filter.Receiver(filtTxSig, spc, chipFilterDelay, BPSes(rate), chipSpreadLengths(rate));
            % demodulate
            RxBits = demodFcns{rate}(RxChips, dataRates(rate));
            % compute number of error bits
            overlapSequenceMsgBin = msgBin(1:end-bitDelay);
            overlapSequenceRxBin = RxBits(bitDelay+1:end);
            % accumulate error and total bits over iterations
            nerror = nerror + sum(overlapSequenceMsgBin ~= overlapSequenceRxBin);
            totalbits = totalbits + (length(RxBits) - bitDelay);
            fprintf('.');
        end
        fprintf('\n');
        BERVector(i, rate) = nerror/totalbits;
    end
end
toc;
fprintf('Simulation completed.\n');
%%
figure;
for rate=1:4
    subplot(2,2,rate);semilogy(snrVector,BERVector(:,rate),'*-');grid on; 
    title_str = {['802.11b: BER for rate = ' num2str(dataRates(rate))...
                 'Mbps'] 'through an AWGN Channel'}; 
    title(title_str); xlabel('Eb/No (dB)'); ylabel('Bit Error Rate');
end
%% 
figure;
for rate=1:4
    semilogy(snrVector, BERVector(:,rate),'^-'); grid on; hold on;
end
hold off;
title('802.11b: BER vs. SNR of Tx through an AWGN channel');
legend('Barker1','Barker2','CCK5.5','CCK11','Location','SouthWest');
xlabel('Eb/No (dB)'); ylabel('Bit Error Rate');