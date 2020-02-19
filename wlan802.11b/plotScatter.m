% plot modulated signal, noisy signal after adding AWGN, and after demod
figure; 
scatter(real(msgMod),imag(msgMod),'filled','m','LineWidth',1.5);
hold on;
scatter(real(txNoisy),imag(txNoisy),'MarkerEdgeColor','r','LineWidth',1);
hold on;
scatter(real(rxBits),imag(rxBits),'MarkerEdgeColor','g','LineWidth',2); hold off;
title('I-Q Phase');
legend('Modulated (Barker1)','After AWGN','Demodulated (Barker1)');