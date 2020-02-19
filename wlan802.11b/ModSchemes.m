classdef ModSchemes
    
    properties(Constant)
        % chipping code length for data rates 1 and 2 Mbps
        %   = static PN (Pseudo Noise) sequence for DSSS
        BarkerSequence = [1 -1  1  1 -1  1  1  1 -1 -1 -1]';
    end
    
    methods(Static)        
        % == DSSS Modulation ==
        function [chippedSignal] = BarkerModulator(data, dataRate)
            % [Usage]
            %   [chippedSignal] = BarkerModulator(data, dataRate)
            % [Parameters]
            %   data: an input sequence of binary bits [0, 1]
            %   dataRate: 1 or 2 Mbps
            % [Output]
            %   chipped code that is spread over Barker code length (11)
            
            % default: initial DBPSK phase: 0.785
            %          initial DQPSK phase: 0
            %          clocksbit: locked
            switch dataRate
                % 1 bit/sec: DBPSK
                case 1
                    modulator = comm.DBPSKModulator;
                    reset(modulator);
                    % modulator.PhaseRotation = 0; % default: 0.785
                % default: ClocksBit: locked
                % 2 bits/sec: DQPSK
                case 2
                    modulator = comm.DQPSKModulator;
                    reset(modulator);
                    modulator.PhaseRotation = 0;
                    modulator.BitInput = true;
                    modulator.SymbolMapping = 'Binary';
                    modulator.OutputDataType = 'double';
            end
            
            % modulate the input data
            % BarkerSequence: [11x1], symbols: [Nx1]
            release(modulator); % clear modulator before stepping
            symbols = modulator(data);
            
            % output: [11 x N/2]
            spreaded = ModSchemes.BarkerSequence*symbols.';
            
            % reshape into a column vector, 
            % each symbol being spreaded over 11 bits
            chippedSignal = spreaded(:);
        end
        function [data] = BarkerDemodulator(chippedSignal, dataRate)
            % [Usage]
            %   [data] = BarkerDemodulator(chippedSignal, dataRate)
            % [Parameters]
            %   chippedSignal: an input sequence of complex chipped codes
            %   dataRate: 1 or 2 Mbps
            % [Output]
            %   binary sequence of data bits, de-spreaded
            
            switch dataRate
                % 1 bit/sec: DBPSK
                case 1
                    % same configuration as modulator
                    modulator = comm.DBPSKDemodulator;
                    reset(modulator);
                    % modulator.PhaseRotation = 0; % default: 0.785
                % 2 bits/sec: DQPSK
                case 2
                    % same configuration as modulator
                    modulator = comm.DQPSKDemodulator;
                    reset(modulator);
                    modulator.PhaseRotation = 0;
                    modulator.BitOutput = true;
                    modulator.SymbolMapping = 'Binary';
                    modulator.OutputDataType = 'double';
            end
            % chippedSignal: (col vector) [(N/2)*11 x 1] --> [11 x (N/2)]
            spreadRate = size(ModSchemes.BarkerSequence,1);
            
            % case for preamble check
            if length(chippedSignal) == 144 || length(chippedSignal) == 48
                appendLength = spreadRate - rem(length(chippedSignal),spreadRate);
                chippedSignal(end+1:end+appendLength) = zeros([appendLength,1]);
            end
            
            chipped = reshape(chippedSignal, spreadRate,[]);
            
            % de-spread the chipped codes
            % symbols: BarkerSequence [1x11] * chipped [11 x (N/2)] = [1 x (N/2)]
            symbols = reshape(ModSchemes.BarkerSequence' * chipped,[],1);
            % make it a column vector, normalize
            symbols = symbols/spreadRate;
            
            % demodulate
            release(modulator);
            data = modulator(symbols);
        end
        
        function [] = plotBarkerSequence(dataRate)
            % plotBarkerSequence(dataRate) plots normalized 
            % autocorrelations of the 11-Barker Sequence 
            % (e.g. power spectrum)
            spreadRate = length(ModSchemes.BarkerSequence);
            acorr = xcorr(ModSchemes.BarkerSequence,'normalized');
            n = linspace(-spreadRate/2-1,spreadRate/2+1,length(acorr));
            figure; plot(n, acorr, 'k-');
            title('11-chip Barker Code Autocorrelation');
            xlabel('Index'); ylabel('normalized autocorr value');
            refline(0,0); legend('Baker Seq 11 chip');
            axis tight; xlim([min(n) max(n)]);
            annotation('textbox',[.15 .3 .3 .15], ...
                'String', strcat('- concentrated main lobe, - very low sidelobes;', ...
                    ' spreading gain:',...
                    num2str(abs(round(10*log10(dataRate/spreadRate),1))),' dB'), ...
                'Color', 'b');
        end
        
        % == High Rate CCK Modulation ==
        % DQPSK encoding of phi_1 from first two bits of data
        % = form of generalized Hadamard encoding
        %
        % * for DQPSK, all odd symbols get an extra 180' rotation
        % Dibit pattern (d0,d1) -- even symbols phase change --   odd sym
        %     00                              0                   pi
        %     01                              pi/2                3pi/2
        %     11                              pi                  0
        %     10                              3pi/2               pi/2
        function [ phi ] = QPSKPhaseTable(dibit)
            % input:  dibit: N rows of data bit pair (column vector)
            %                assume EVEN symbols
            % output: column vector of N values for each row dibit input

            % returns one phi value for the input pair, for each row
            % phase assignment per phase angle PHI
            phi = zeros(size(dibit,1), 1);
            for i = 1:length(phi)
                switch char(dibit(i,:))
                    case char([0 0])
                        phi(i) = 0;
                    case char([0 1])
                        phi(i) = pi/2;
                    case char([1 0])
                        phi(i) = pi;
                    case char([1 1])
                        phi(i) = 3*pi/2;
                end
            end
        end
        
        function [ phi ] = DQPSK_PHI1(dibit)
            % input: dibit for first phase angle -- PHI 1
            % output: offsetted phases for all N rows (DQPSK)
            
            % pi offset vector
            pioffset = zeros([size(dibit,1),1]);
            for jj = 1:size(dibit,1)
                if mod(jj,2) == 1       % even
                    pioffset(jj) = 0;
                elseif mod(jj,2) == 0   % odd
                    pioffset(jj) = pi;
                end
            end
                  
            % add to previous symbol phase angle
            phi = cumsum(ModSchemes.QPSKPhaseTable(dibit)) + pioffset;
            
        end
            
        function [phi] = CCKPhases5_5(data_4bit)
            % output: each phi_i is a column vector; 
            %         phi is a matrix of 4 phi's, each in its column
            % phi 1 is "even" and assume previous phase angle is 0
            phi_1 = ModSchemes.DQPSK_PHI1(data_4bit(:,1:2));
            phi_2 = data_4bit(:,3) * pi + pi/2;
            phi_3 = zeros(size(phi_1));
            phi_4 = data_4bit(:,4) * pi;
            phi = [phi_1, phi_2, phi_3, phi_4];
        end
        
        function [phi] = CCKPhases11(data_8bit)
            % output: each phi_i is a column vector; 
            %         phi is a matrix of 4 phi's, each in its column
            phi_1 = ModSchemes.DQPSK_PHI1(data_8bit(:,1:2));
            phi_2 = ModSchemes.QPSKPhaseTable(data_8bit(:,3:4));
            phi_3 = ModSchemes.QPSKPhaseTable(data_8bit(:,5:6));
            phi_4 = ModSchemes.QPSKPhaseTable(data_8bit(:,7:8));
            phi = [phi_1, phi_2, phi_3, phi_4];
        end
        
        function [CCKsymbols] = CCKWordGenerate(phi)
            % outputs the codeword for the specified phase angles PHI 1-4
            phi_1 = phi(:,1);
            phi_2 = phi(:,2);
            phi_3 = phi(:,3);
            phi_4 = phi(:,4);
            CCKsymbols = [exp(1j*(phi_1+phi_2+phi_3+phi_4)), ...
                         exp(1j*(phi_1+phi_3+phi_4)), ... 
                         exp(1j*(phi_1+phi_2+phi_4)), ...
                         -exp(1j*(phi_1+phi_4)), ... % 180' offset
                         exp(1j*(phi_1+phi_2+phi_3)), ...
                         exp(1j*(phi_1+phi_3)), ...
                         -exp(1j*(phi_1+phi_2)), ... % 180' offset
                         exp(1j*(phi_1))
                        ];
        end
    
        function [symbols] = CCKModulator(data, dataRate)
            % input: data bit stream (a single vector) and data rate
            % output: modulated symbol stream for the specified data rate
            
            % set bps to truncate by
            if nargin < 2
                % default: 11 Mbps / 8 bits per symbol
                dataRate = 11;
                bps = 8; % bits per symbol
            else
                % if dataRate is specified, set bps
                if dataRate == 5.5
                    bps = 4;
                elseif dataRate == 11
                    bps = 8;
                end
            end
            
            % add zeros to make it multiple of bps
            if mod(length(data),bps) > 0
                data = [data(:); zeros([bps - mod(length(data),bps),1])];
            end
            data = reshape(data,bps,[])';
            
            % perform modulation by data rate
            switch dataRate
                case 5.5
                    % 4 bits per sample
                    phi = ModSchemes.CCKPhases5_5(data);
                case 11
                    % 8 bits per sample
                    phi = ModSchemes.CCKPhases11(data);
            end
            % generate a length-8 chipping codeword from 4 phase angles
            complexChips = ModSchemes.CCKWordGenerate(phi).';
            symbols = complexChips(:); % make a column vector
%             scatterplot(complexChips(:))
        end
        
        function [rxBits] = CCKDemodulator(inputSymbols, dataRate)
            % takes in symbol stream (chipped) and demodulates to obtain a
            % bit data stream in binary
            % * bitDelay is set to be a multiple of bps
            
            if dataRate == 5.5
                bps = 4;
            elseif dataRate == 11
                bps = 8;
            end

            symsOctet = reshape(inputSymbols, 8, []).';

            % generate CCK modulated symbols for all states - for phase angles
            codewordsPHI23 = zeros(8, 2^(bps-2));
            % set first dibit to a random init value for 8-bit data input 
            % for CCK symbols below, to be used to search for states
            initDibit = [0, 0]; 
            if dataRate == 11
                % data rate of 11 Mbps / 8 bits per symbol
                for n = 1:2^6
                    bits2to7 = de2bi(n-1, bps-2, 'left-msb');
                    % returns each codeword in a column vector,
                    % 64 of them in each column
                    % 64: total number of possible phase 2-4 states to search for phi 1
                    codewordsPHI23(:,n) = ModSchemes.CCKModulator([initDibit, bits2to7], 11);
                end
            elseif dataRate == 5.5
                % data rate of 5.5 Mbps / 4 bits per symbol
                for n = 1:2^2
                    bits2to7 = de2bi(n-1, bps-2, 'left-msb');
                    % returns each codeword in a column vector,
                    % 4 of them in each column
                    % 4: total number of possible phase 2-4 states to search for phi 1
                    codewordsPHI23(:,n) = ModSchemes.CCKModulator([initDibit, bits2to7], 5.5);
                end
            end

            % (binary) table of possible states for dibit corresponding to PHI 1
            dibitStates = [0, 0; 
                           0, 1; 
                           1, 0; 
                           1, 1]; %possibilities for first 2 bits

            % initial prediction for PHI 1
            % note that 5.5 Mbps data rate is just QPSK, while
            % 11 Mbps uses DQPSK, which depends on its prior value
            estPhi1 = 1;
                       
            % separate complex symbols into real and imaginary components
            allcodewordsPHI23 = [real(codewordsPHI23); imag(codewordsPHI23)];

            estIdx = zeros([size(symsOctet,1), 1]);
            for n = 1:size(symsOctet, 1) % loop through each of N rows
                % c7 = phi 1 estimate
                % obtain real phase angle for PHI 1 by multiplying
                % conjugates of complex codewords with another complex
                % codeword (c7)
                c7_PHI1_codeword_est = symsOctet(n,:) * conj(symsOctet(n, 8)); 
                realimag_c7_PHI1_codeword_est = [real(c7_PHI1_codeword_est), imag(c7_PHI1_codeword_est)];

                % for odd symbols, add extra 180' shift
                % estPhi1 is running estimate of phi 1 (based on previous values)
                evenOddIndex = n - 1;
                if (mod(evenOddIndex,2)) % even
                    % no shift if even symbols --> 0'
                    noshift = exp(0);
                    shiftedPHI1 = symsOctet(n, 8) * conj(estPhi1) * noshift;
                else % odd
                    % shift for DQPSK odd symbols (shift == 1) --> 180'
                    pishift = exp(-1*pi*1j);
                    shiftedPHI1 = symsOctet(n, 8) * conj(estPhi1) * pishift;
                end
                
                % every chip in a code word can receive 4 values {1,-1, j,-j}, 
                % so there are a total of 4^8 = 65536 different combinations of 8 chips.
                allValuesForChip = [1, 1j, -1, -1j]; % in a codeword c0-07
                
                % find the shifted PHI 1 value from the table of 4 states
                % dsearchn: N-dimensional nearest point search for PHI 1 
                %            without using a triangulation
                %           - returns the indices fo the closest points
                idx = dsearchn([real(allValuesForChip); imag(allValuesForChip)].', ...
                                [real(shiftedPHI1), imag(shiftedPHI1)]);
                dibit_1 = dibitStates(idx, :);
                % dibit represented in 8-bit
                cumulative = bi2de([dibit_1, zeros([1, bps-2])], 'left-msb');
                
                % add the odd (pi) shift to the even phase angle
                idxoffset = 1; % to start from 0
                estIdx(n,1) = dsearchn(allcodewordsPHI23.', realimag_c7_PHI1_codeword_est) + cumulative - idxoffset;
                
                estPhi1 = symsOctet(n, 8); % update estimate with c7 value
            end
            % make column vector
            rxBits = reshape(de2bi(estIdx, bps, 'left-msb').', [], 1); 
        end
    end % end methods
end % end classdef
