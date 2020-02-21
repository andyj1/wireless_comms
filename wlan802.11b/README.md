##### WLAN IEEE 802.11b

- Simulation of IEEE 802.11b (Wireless Local Area Network) in MATLAB 
  - (R2019b version; Toolbox: communications, DSP (for modulator), WLAN (for scrambler)
- For code, _src/wlan80211b_packetframe.m_ should contain packet frame transmission (preamble, header verification between Tx and Rx is performed only)
- Possible improvements: better CCK decoding technique, RAKE receiver with equalizers, etc.

##### IEEE 802.11b

- Modulation schemes: DSSS (1, 2 Mbps), CCK (5.5, 11 Mbps)
- Coding sequences: Barker-11 chipping sequence (static PN sequence), CCK codeword
- Pulse Shaping Filter: Root Raised Cosine Fitler ($M$ = 40, $\beta$ = 0.3)

> ##### Glossary
>
> - DSSS: Direct Sequence Spread Spectrum
> - CCK: Complementary Code Keying

##### Class definitions

- _Filter_: Tx and Rx along with up/down-sampling
- _ModSchemes_: BarkerModulator/Demodulator, CCKModulator/Demodulator, plotBarkerAutocorrelation

> A written report is available [here](doc/ece408_802.11b_report.pdf)

> A published file (MATLAB) is available [here](doc/wlan80211b_packetframe.pdf)

##### Bit Error Rate Curves

![BER curve 1](res/ber_all_rates.png)
![BER curve 2](res/ber_per_rate.png)
