if ~mod(tt,plot_int)
    figure(1)
    subplot 321
    plot(Grid.xc,S_new(:,2),'r-','Displayname','Saq','linewidth',3)
    hold on
    plot(Grid.xc,S_new(:,3),'g-','Displayname','Sv','linewidth',3)
    plot(Grid.xc,S_new(:,1),'b-','Displayname','Sh','linewidth',3)
    xlabel('Distance (x)')
    ylabel('Saturation')
    legend('show')
    xlim([0 1])
    ylim([0 1])
    title(['time =',num2str(tt*dt)])
    hold off

    subplot 322
    aq_mole = squeeze(x_new(:,1,:));
    aq_gst_mole = aq_mole(2:end,:)./repmat(sum(aq_mole(2:end,:),1), [3 1])./repmat((S_new(:,1)>0)', [3 1]);
    plot(Grid.xc,aq_gst_mole(1,:),'r-','Displayname','CH4','linewidth',3)
    hold on
    plot(Grid.xc,aq_gst_mole(2,:),'b-','Displayname','CO2','linewidth',3)
    plot(Grid.xc,aq_gst_mole(3,:),'g-','Displayname','N2','linewidth',3)
    xlabel('Distance (x)')
    ylabel('Aqueous guest mole fraciotn')
    legend('show')
    xlim([0 1])
    ylim([0 1])
    title(['time =',num2str(tt*dt)])
    hold off

    subplot 323
    v_mole = squeeze(x_new(:,2,:));
    v_gst_mole = v_mole(2:end,:)./repmat(sum(v_mole(2:end,:),1), [3 1])./repmat((S_new(:,2)>0)', [3 1]);
    plot(Grid.xc,v_gst_mole(1,:),'r-','Displayname','CH4','linewidth',3)
    hold on
    plot(Grid.xc,v_gst_mole(2,:),'b-','Displayname','CO2','linewidth',3)
    plot(Grid.xc,v_gst_mole(3,:),'g-','Displayname','N2','linewidth',3)
    xlabel('Distance (x)')
    ylabel('Vapor guest mole fraciotn')
    legend('show')
    xlim([0 1])
    ylim([0 1])
    title(['time =',num2str(tt*dt)])
    hold off

    subplot 324
    h_mole = squeeze(x_new(:,3,:));     
    h_gst_mole = h_mole(2:end,:)./repmat(sum(h_mole(2:end,:),1), [3 1])./repmat((S_new(:,3)>0)', [3 1]);
    plot(Grid.xc,h_gst_mole(1,:),'r-','Displayname','CH4','linewidth',3)
    hold on
    plot(Grid.xc,h_gst_mole(2,:),'b-','Displayname','CO2','linewidth',3)
    plot(Grid.xc,h_gst_mole(3,:),'g-','Displayname','N2','linewidth',3)
    xlabel('Distance (x)')
    xlabel('Distance (x)')
    ylabel('Hydrate cage percent occupied')
    legend('show')
    xlim([0 1])
    ylim([0 1])
    title(['time =',num2str(tt*dt)])
    hold off

    subplot 325
    plot(Grid.xc,vd,'k-','linewidth',3)
    xlabel('Distance (x)')
    ylabel('Local Velocity')
    ylim([0 3])
    xlim([0 1])
    title(['time =',num2str(tt*dt)])
    hold off

    subplot 326
    plot(Grid.xc,rho,'k-','DisplayName','rho','linewidth',3)
    hold on
    plot(Grid.xc,sum(G_new,2),'b--','Displayname','G','linewidth',3)
    xlabel('Distance (x)')
    ylabel('Density')
    xlim([0 1])
    title(['time =',num2str(tt*dt)])
    legend('show')
    drawnow
    hold off
end