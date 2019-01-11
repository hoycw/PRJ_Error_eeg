function fn_plot_PSD_1by1_save(data1, labels, sample_freq, save_prefix, save_filetype)
% Check noise profile
for channel_n = 1:length(labels)
    figure('Visible','off');
    [fft_data,freqs] = pwelch(data1(channel_n,:),2048,0,2048,sample_freq);
    loglog(freqs,fft_data);
    xlim([1 350]);
    ax = gca;
    ax.XTick = [4 8 12 25 30 60 80 100 120 180 200 240 300 360 420];
    title(['Channel ' num2str(channel_n) ' [' labels{channel_n} ']']);
    saveas(gcf,strcat(save_prefix,'_',labels{channel_n},'.',save_filetype));
    close;
end
end
