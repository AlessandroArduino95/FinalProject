%% PLOTS USED

plot(t_val(:),u_time((length(nodes)/2 + 0.5),:),'-','LineWidth',8) % I take the values from u(2) forward as u(1) is the stationary solution
title('Evolution of the solution over time S=1')
xlabel('time')
ylabel('u (solution)')
set(gca,'fontsize', 27);


%initial
u_toplot = reshape(u_time(:,1),length(x_val),[]);
% length(t_val)
figure
surf(x_val,y_val,u_toplot,'EdgeColor','interp')
title('Initial solution (5x5 grid)')
xlabel('x')
ylabel('y')
set(gca,'fontsize', 33)

%final
u_toplot = reshape(u_time(:,length(t_val)),length(x_val),[]);
% length(t_val)
figure
surf(x_val,y_val,u_toplot,'EdgeColor','interp')
title('Final solution (5x5 grid)')
xlabel('x')
ylabel('y')
set(gca,'fontsize', 33);
%zlim([0 5]);