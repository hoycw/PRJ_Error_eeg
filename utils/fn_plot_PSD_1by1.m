function fn_plot_PSD_1by1(data, labels, sample_freq)
% Check noise profile
for channel_n = 1:length(labels)
    [fft_data,freqs] = pwelch(data(channel_n,:),2048,0,2048,sample_freq);
    loglog(freqs,fft_data);
    xlim([1 350]);
    ax = gca;
    ax.XTick = [4 8 12 25 30 60 80 100 120 180 200 240 300 360 420];
    title(['Channel ' num2str(channel_n) ' [' labels{channel_n} ']']);
    pause;
end
end
