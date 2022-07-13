function coherence_plots(brain_out,aal_out,plot_folder,freqs)

figure
set(gcf,'position',[1 1 1440 700])
tiledlayout(2,5,'tilespacing','tight','padding','tight')
fig_name = 'Fig S2';
fid = fopen([plot_folder,'supplemental_results.html'],'a');
fprintf(fid,'<br><b>Coherence-based connectivity with symmetric coverage constraint</b></br>');

freq_letters = {'A','B','C','D','E','F','G','H','I','J'};

for ia = 1:2
    
    if ia == 1
        atlas_out = brain_out;
        fprintf(fid,'<p><br><b>Brainnetome atlas</b></br>');
        fprintf(fid,['First, using the Brainnetome atlas to parcellate brain regions,',...
            ' we compared the frequency-specific coherence to the SOZ versus '...
            'that to the contralateral region.']);
        
    else
        atlas_out = aal_out;  
        fprintf(fid,'<br><b>AAL atlas</b></br>');
        fprintf(fid,['Next, using the AAL atlas to parcellate brain regions,',...
            ' we compared the frequency-specific coherence to the SOZ versus '...
            'that to the contralateral region.']);
    end
    
    % unpack brain_out - be careful not to confuse with aal results
    unpack_any_struct(atlas_out);
    
    for f = 1:5
        
        nexttile
        if ia == 1 && f == 1
            stats = paired_plot(squeeze(soz_coh_all(:,f,:)),'Coherence',...
                {'SOZ','in contralateral region','in contralateral region'},...
                0,0,[0.045,0.80,0.16,0.08]);

        else
            stats = paired_plot(squeeze(soz_coh_all(:,f,:)),...
                'Coherence',{'SOZ','in contralateral region','in contralateral region'},...
                0,1);
        end
        if ia == 1
            title(sprintf('Brainnetome\n%s',freqs{f}))
        else
            title(sprintf('AAL\n%s',freqs{f}))
        end
        if stats.pval < 0.05
            fprintf(fid,[' The average %s coherence to the SOZ (median %s) was lower than that to '...
                ' the region contralateral to the SOZ (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%s%s).'],...
                freqs{f},pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name,freq_letters{f+(ia-1)*5});
        else
            fprintf(fid,[' There was no significant difference between the %s coherence to the SOZ (median %s) and that to '...
                ' the region contralateral to the SOZ (median %s) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (%s%s).'],...
                freqs{f},pretty_exp_html(stats.medians(1)),pretty_exp_html(stats.medians(2)),stats.Tpos,get_p_html(stats.pval),fig_name,freq_letters{f+(ia-1)*5});
        end
        
    end
    
    fprintf(fid,'</p>');

end

% Add annotations
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.2 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0.4 0.91 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.6 0.91 0.1 0.1],'String','D','fontsize',25,'linestyle','none')
annotation('textbox',[0.8 0.91 0.1 0.1],'String','E','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.41 0.1 0.1],'String','F','fontsize',25,'linestyle','none')
annotation('textbox',[0.2 0.41 0.1 0.1],'String','G','fontsize',25,'linestyle','none')
annotation('textbox',[0.4 0.41 0.1 0.1],'String','H','fontsize',25,'linestyle','none')
annotation('textbox',[0.6 0.41 0.1 0.1],'String','I','fontsize',25,'linestyle','none')
annotation('textbox',[0.8 0.41 0.1 0.1],'String','J','fontsize',25,'linestyle','none')
print(gcf,[plot_folder,fig_name],'-dpng')
fclose(fid);

end