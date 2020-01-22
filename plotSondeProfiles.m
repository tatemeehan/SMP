        % plot Rut Profiles
        figure();
        subplot(1,3,1)
        for kk = 1:length(smprutix)
            if kk == 1
                smpmax = max(invSMP{smprutix(kk)}.M(:,7));
            elseif max(invSMP{smprutix(kk)}.M(:,7)) > smpmax
                smpmax = max(invSMP{smprutix(kk)}.M(:,7));
            end
        end
                
        for kk = 1:length(smprutix)
            plot(invSMP{smprutix(kk)}.M(:,7)./smpmax,invSMP{smprutix(kk)}.z,'r','linewidth',1.5); hold on;
        end
        rammbut = rammrut(rammrut ~= 46);
        subplot(1,3,1)
        for kk = 1:length(rammbut)
            if kk == 1
                rammmax = max(ramm(rammbut(kk)).index(:));
            elseif max(ramm(rammbut(kk)).index(:)) > rammmax
                rammmax = max(ramm(rammbut(kk)).index(:));
            end
        end
        for kk = 1:length(rammbut)
            stairs(ramm(rammbut(kk)).index(:)./rammmax,ramm(rammbut(kk)).penetration*10,'k','linewidth',2);hold on
%             plot(ramm(rammbut(kk)).index(:)./rammmax,ramm(rammbut(kk)).penetration*10,'k');hold on
        end
        ylabel('Depth [mm]')
        set(gca,'fontweight','bold','fontsize',14)
                title({' ','Rut'})
        axis ij
        % Plot Belly Profiles
        subplot(1,3,2)
        for kk = 1:length(smpbellyix)
            if kk == 1
                smpmax = max(invSMP{smpbellyix(kk)}.M(:,7));
            elseif max(invSMP{smpbellyix(kk)}.M(:,7)) > smpmax
                smpmax = max(invSMP{smpbellyix(kk)}.M(:,7));
            end
        end
        for kk = 1:length(smpbellyix)
            if kk == length(smpbellyix)
            hkk(1) = plot(invSMP{smpbellyix(kk)}.M(:,7)./smpmax,invSMP{smpbellyix(kk)}.z,'r','linewidth',1.5); hold on;
            else
                plot(invSMP{smpbellyix(kk)}.M(:,7)./smpmax,invSMP{smpbellyix(kk)}.z,'r','linewidth',1.5); hold on;
            end
        end
        for kk = 1:length(rammbelly)
            if kk == 1
                rammmax = max(ramm(rammbelly(kk)).index(:));
            elseif max(ramm(rammbelly(kk)).index(:)) > rammmax
                rammmax = max(ramm(rammbelly(kk)).index(:));
            end
        end
        for kk = 1:length(rammbelly)
            if kk == length(rammbelly)
            hkk(2) = stairs(ramm(rammbelly(kk)).index(:)./rammmax,ramm(rammbelly(kk)).penetration*10,'k','linewidth',2);hold on
            else
            hkk(2) = stairs(ramm(rammbelly(kk)).index(:)./rammmax,ramm(rammbelly(kk)).penetration*10,'k','linewidth',2);hold on
            end                
%             plot(ramm(rammbelly(kk)).index(:)./rammmax,ramm(rammbelly(kk)).penetration*10,'k');hold on
        end
        legend([hkk],'SMP','RammSonde','location','southeast')
        xlabel('Relative Strength (Index)')
        set(gca,'fontweight','bold','fontsize',14)
        title({SVtag{1}{8},'Belly'})
        axis ij
        % Plot VS Profiles
        for kk = 1:length(smpvsix)
            if kk == 1
                smpmax = max(invSMP{smpvsix(kk)}.M(:,7));
            elseif max(invSMP{smpvsix(kk)}.M(:,7)) > smpmax
                smpmax = max(invSMP{smpvsix(kk)}.M(:,7));
            end
        end        
         subplot(1,3,3)
        for kk = 1:length(smpvsix)
        plot((invSMP{smpvsix(kk)}.M(:,7)./smpmax),invSMP{smpvsix(kk)}.z,'r','linewidth',1.5); hold on;
        ylim([0,550])
        end
        for kk = 1:length(rammvs)
            if kk == 1
                rammmax = max(ramm(rammvs(kk)).index(:));
            elseif max(ramm(rammvs(kk)).index(:)) > rammmax
                rammmax = max(ramm(rammvs(kk)).index(:));
            end
        end
        for kk = 1:length(rammvs)
            stairs(ramm(rammvs(kk)).index(:)./rammmax,ramm(rammvs(kk)).penetration*10,'k','linewidth',2);hold on
%             plot(ramm(rammvs(kk)).index(:)./rammmax,ramm(rammvs(kk)).penetration*10,'k');hold on
        end
        set(gca,'fontweight','bold','fontsize',14)
                title({' ','Virgin Snow'})
        axis ij