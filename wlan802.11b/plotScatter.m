function plotScatter(fig, msgMod, symbolsTx, rxSyms, rate)
    % plot modulated signal, noisy signal after adding AWGN, and after demod
    
    scatter(real(symbolsTx),imag(symbolsTx),50,'filled','y','MarkerEdgeColor','y','LineWidth',3);
    hold on;
    scatter(real(msgMod),imag(msgMod),50,'filled','r');
    hold on;
    scatter(real(rxSyms),imag(rxSyms),5, 'filled','g'); 
    hold off;
    title(strcat(num2str(rate),' Mbps Scatterplot for PPDU'));
    xlabel('In-phase amplitude'); ylabel('Quadrature Amplitude');
    xlim([-1.1 1.1]); ylim([-1.1 1.1]);
    legend('Tx Symbols','AWGN Channel','Rx Symbols');
    set(gca, 'FontSize', 12, 'FontWeight','bold');  
end