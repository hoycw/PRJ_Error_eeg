function fn_plot_PSD_1by1_compare_save(data1, data2, labels1, labels2, sample_freq,...
                                       save_prefix, data_type1, data_type2, save_filetype)

if numel(labels1)~=numel(labels2)
    fprintf('!!!WARNING!!! different number of channels in data1 (%i) and data2 (%i)\n',...
        numel(labels1),numel(labels2));
end
min_n_labels = min(numel(labels1),numel(labels2));

% Check noise profile
for channel_n = 1:min_n_labels
    figure('Visible','off');
    [fft_data1,freqs1] = pwelch(data1(channel_n,:),2048,0,2048,sample_freq);
    [fft_data2,freqs2] = pwelch(data2(channel_n,:),2048,0,2048,sample_freq);
    loglog(freqs1,fft_data1,'k');
    hold on;
    loglog(freqs2,fft_data2,'r');
    xlim([1 350]);
    ax = gca;
    ax.XTick = [4 8 12 25 30 60 80 100 120 180 200 240 300 360 420];
    title(['Channel ' num2str(channel_n) ' [' labels1{channel_n} 'vs' labels2{channel_n} ']']);
    legend(strcat(data_type1,'-',labels1{channel_n}), strcat(data_type2,'-',labels2{channel_n}));
    saveas(gcf,strcat(save_prefix,'_',labels1{channel_n},'.',labels2{channel_n},'.',save_filetype));
    close;
end
end
