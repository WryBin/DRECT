function contourDOSY(Type, Result_LRsPILT, Result_CoMeF, Z, HNMR, b, ppm, cs, idx_peaks, cs1, cs2, dc1, dc2, t)
    cs_spec = zeros([length(ppm), 1]);
    spec_whole = zeros([length(Z(1, :)), length(ppm)]);
    cs_spec(idx_peaks, :) = HNMR;
    spec_whole(:, idx_peaks) = Z.';
    decay_range = linspace(0, (length(Z(1, :))-1)/10, length(Z(1, :)));

    if Type == "VD"
        load VD_NMR_spectra.mat
        plot(ppm,NMR_spectra, "Color",'k');set(gca,'Xdir','reverse');axis off;
        xlim([cs1,cs2]);

        % LRsPLIT
        nexttile(4, [2, 1])
        DiffCoef = [13.23, 13.92];
        contour(cs,linspace(1, 20, 191), Result_LRsPILT, 40);

        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([dc1,dc2]);
        text(2.47, 12.7, 'b', FontSize=12)
        text(2.38, 13.23, "vitamin D_3", FontSize=8, Color=[1, 0.07, 0.07])
        text(2.38, 13.92, "provitamin D_3", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );
        set(gca,'xtick',[])

        
        % CoMeF
        nexttile(6, [2, 1])
        DiffCoef = [12.3, 12.7];
        contour(linspace(1.8015, 2.4996, 287),linspace(12, 14, 100), Result_CoMeF, 20);

        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([12, 13]);
        text(2.47, 12.1, 'c', FontSize=12)
        text(2.38, 12.3, "vitamin D_3", FontSize=8, Color=[1, 0.07, 0.07])
        text(2.38, 12.7, "provitamin D_3", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );
        set(gca,'xtick',[])

        
        % DRILT
        nexttile(8, [2, 1])
        DiffCoef = [13.23, 13.92];
        contour(ppm,(decay_range)/7 + 13,spec_whole,40);
        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([12.7,14.3]);
        text(2.47, 13.0, 'd', FontSize=12)
        text(2.38, 13.23, "vitamin D_3", FontSize=8, Color=[1, 0.07, 0.07])
        text(2.38, 13.92, "provitamin D_3", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );

        xlabel(t, 'Chemical Shift(ppm)');
        ylabel(t, 'Diffusion Coefficient(10^{-10}m^2/s)');

    elseif Type == "GSP"
        plot(ppm,cs_spec, "Color",'k');set(gca,'Xdir','reverse');axis off;
        xlim([cs1,cs2]);

        % LRsPLIT
        nexttile(2, [2, 1])
        DiffCoef = [2.2, 3.12, 4.02];
        contour(ppm,linspace(0, 7, 351), Result_LRsPILT, 20);

        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([dc1,dc2]);
        text(5.4, 1.9, 'a', FontSize=12)
        text(5.0, 2.2, "PEG600", FontSize=8, Color=[1, 0.07, 0.07])
        text(5.0, 3.07, "sucrose", FontSize=8, Color=[1, 0.07, 0.07])
        text(5.0, 3.97, "glucose", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );
        set(gca,'xtick',[])

        
        % CoMeF
        nexttile(4, [2, 1])
        DiffCoef = [2.2, 3.11, 4.01];
        contour(ppm,linspace(1.5, 5, 100), Result_CoMeF, 20);

        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([dc1,dc2]);
        text(5.4, 1.9, 'b', FontSize=12)
        text(5.0, 2.2, "PEG600", FontSize=8, Color=[1, 0.07, 0.07])
        text(5.0, 3.07, "sucrose", FontSize=8, Color=[1, 0.07, 0.07])
        text(5.0, 3.97, "glucose", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );
        set(gca,'xtick',[])

        
        % DRILT
        nexttile(6, [2, 1])
        DiffCoef = [2.21, 3.02, 4.14];
        contour(ppm,decay_range*(0.8/b(end)),spec_whole, 40);
        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([dc1,dc2]);
        text(5.4, 1.9, 'c', FontSize=12)
        text(5.0, 2.16, "PEG600", FontSize=8, Color=[1, 0.07, 0.07])
        text(5.0, 2.97, "sucrose", FontSize=8, Color=[1, 0.07, 0.07])
        text(5.0, 4.09, "glucose", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );

        xlabel(t, 'Chemical Shift(ppm)');
        ylabel(t, 'Diffusion Coefficient(10^{-10}m^2/s)');

    elseif Type == "QGC"
        plot(ppm,cs_spec, "Color",'k');set(gca,'Xdir','reverse');axis off;
        xlim([cs1,cs2]);

        % LRsPILT
        nexttile(2, [2, 1])
        DiffCoef = [4.7, 7.3, 10.1];
        contour(ppm,linspace(1, 20, 191), Result_LRsPILT, 40);

        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([dc1,dc2]);
        text(12.2, 3, 'a', FontSize=12)
        text(10.1, 4.66, "quinine", FontSize=8, Color=[1, 0.07, 0.07])
        text(10.1, 7.26, "geraniol", FontSize=8, Color=[1, 0.07, 0.07])
        text(10.1, 10.06, "camphene", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );
        set(gca,'xtick',[])

        % CoMeF
        nexttile(4, [2, 1])
        DiffCoef = [4.6, 7.3, 9.9];
        contour(ppm,linspace(1, 12, 100), Result_CoMeF, 40);

        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([dc1,dc2]);
        text(12.2, 3, 'b', FontSize=12)
        text(10.1, 4.56, "quinine", FontSize=8, Color=[1, 0.07, 0.07])
        text(10.1, 7.26, "geraniol", FontSize=8, Color=[1, 0.07, 0.07])
        text(10.1, 9.86, "camphene", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );
        set(gca,'xtick',[])

        % DRILT
        nexttile(6, [2, 1])
        DiffCoef = [4.6, 7.2, 10.4];
        contour(ppm,decay_range*(0.8/b(end)),spec_whole,40);

        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([dc1,dc2]);
        text(12.2, 3, 'c', FontSize=12)
        text(10.1, 4.46, "quinine", FontSize=8, Color=[1, 0.07, 0.07])
        text(10.1, 7.16, "geraniol", FontSize=8, Color=[1, 0.07, 0.07])
        text(10.1, 10.36, "camphene", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );

        xlabel(t, 'Chemical Shift(ppm)');
        ylabel(t, 'Diffusion Coefficient(10^{-10}m^2/s)');

    elseif Type == "M6"
        plot(ppm,cs_spec, "Color",'k');set(gca,'Xdir','reverse');axis off;
        xlim([cs1,cs2]);

        
        % LRsPILT
        nexttile(2, [2, 1])
        DiffCoef = [3.2, 4.2, 5.2, 6.3, 8.2, 10.6];
        contour(ppm,linspace(1, 20, 191), Result_LRsPILT, 40);

        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([dc1,dc2]);
        text(5.3, 3, 'a', FontSize=12)
        text(5, 3.16, "sucrose", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 4.16, "lysine", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 5.16, "threonine", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 6.26, "sucrose", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 8.16, "lysine", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 10.56, "threonine", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );
        set(gca,'xtick',[])

        % CoMeF
        nexttile(4, [2, 1])
        DiffCoef = [3.1, 3.2, 3.3, 4.2, 4.3, 4.9, 6.2, 8.2, 10.7];
        contour(cs,linspace(2, 12, 100), Result_CoMeF, 20);

        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([dc1,dc2]);
        text(5.3, 3, 'b', FontSize=12)
        text(5, 3.1, "sucrose", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 4.16, "lysine", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 4.86, "threonine", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 6.16, "sucrose", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 8.16, "lysine", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 10.66, "threonine", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );
        set(gca,'xtick',[])

        % DRILT
        nexttile(6, [2, 1])
        DiffCoef = [3.2, 4.3, 5.0, 6.3, 8.2, 10.5];
        contour(ppm,decay_range*(1.2/b(end)),spec_whole,40);

        set(gca,'Ydir','reverse','Xdir','reverse'); 
        xlim([cs1,cs2]);
        ylim([dc1,dc2]);
        text(5.3, 3, 'c', FontSize=12)
        text(5, 3.16, "sucrose", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 4.26, "lysine", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 4.96, "threonine", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 6.26, "sucrose", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 8.16, "lysine", FontSize=8, Color=[1, 0.07, 0.07])
        text(5, 10.46, "threonine", FontSize=8, Color=[1, 0.07, 0.07])

        for i = 1:length(DiffCoef)
            l = line(gca,get(gca,'xlim'),DiffCoef(i)*ones(1,2),'LineWidth',0.8,'color',[0.85 0.85 0.85],'LineStyle','--');
            uistack(l, "bottom");
        end
        set(gca,'YTick',unique(DiffCoef) );

        xlabel(t, 'Chemical Shift(ppm)');
        ylabel(t, 'Diffusion Coefficient(10^{-10}m^2/s)');
    end

end