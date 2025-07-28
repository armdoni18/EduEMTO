function [] = Func9_Opt_Plot(fem,opt)

%% 1. Objective function history
    figure(1); clf(1); h=figure(1);
    set(h, 'Position', [1020, 50, 400, 350]);
    plot(opt.fhis,'-k','linewidth',2);
    grid on;
    xlabel('Iteration');
    ylabel('Objective');

%% 2. Volume Constraint function history
    figure(2); clf(2); h=figure(2);
    set(h, 'Position', [1425, 50, 400, 350]);
    plot(opt.ghis + opt.volfrac,'-b','linewidth',2);
    hold on;
    xlabel('Iteration');
    ylabel('Volume');
    plot([1 opt.iter],[opt.volfrac opt.volfrac],':b','linewidth',2); grid on; ylim([0 1]);

%% 3. Element-wise Projected Density
    figure(3); clf(3); h=figure(3);
    set(h, 'Position', [1220, 485, 600, 350]);
    low_density_mask = opt.erho(opt.erho_dd) > 0.99;
    patch('Faces', fem.IX(low_density_mask,1:3), ...
          'Vertices', fem.X, ...
          'FaceVertexCData', opt.erho(low_density_mask), ...
          'FaceColor', 'flat', ...
          'EdgeColor', 'none');
    colorbar;
    xlabel('X'); ylabel('Y');
    title('Element-wise Projected Density (\rho_{element})');
    colormap(flipud(gray));
    axis equal;
    xlim([min(fem.X(:,1)), max(fem.X(:,1))]);
    ylim([min(fem.X(:,2)), max(fem.X(:,2))]);
    saveas(h, sprintf('Figures/Density_Iter_%03d.png', opt.iter));

%% 4. Vector Potential Contours
    x_nodes = fem.X(:,1);                   
    y_nodes = fem.X(:,2);                   
    A_values = fem.A;                       

    [Xq, Yq] = meshgrid(linspace(min(x_nodes), max(x_nodes), 391), ...
                         linspace(min(y_nodes), max(y_nodes), 251));
    Aq = griddata(x_nodes, y_nodes, A_values, Xq, Yq, 'cubic');

    figure(4); clf(4); h=figure(4);
    set(h, 'Position', [1830, 50, 600, 350]);
    contour(Xq, Yq, Aq,15, "k");
    axis equal;
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Vector Potential A_z Contours');
    grid on; drawnow
    saveas(h, sprintf('Figures/A_Iter_%03d.png', opt.iter));

%% 5. Magnetic Flux Density |B|
    figure(5); clf(5); h=figure(5);
    set(h, 'Position', [1830, 485, 600, 350]);
    patch('Faces', fem.IX(:,1:3), ...
          'Vertices', fem.X, ...
          'FaceVertexCData', fem.B, ...
          'FaceColor', 'flat', ...
          'EdgeColor', 'none');
    axis equal; colorbar;
    xlabel('X'); ylabel('Y');
    title('Magnetic Flux Density |B|');
    colormap(parula);
    saveas(h, sprintf('Figures/B_Iter_%03d.png', opt.iter));

    % Optional: print time every 10 iterations
    if rem(opt.iter,10)==0; toc; end    

end
