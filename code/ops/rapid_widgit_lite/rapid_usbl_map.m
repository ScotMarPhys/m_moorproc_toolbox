function rapid_usbl_map(web_dis,work_dir,this_moor,mlat,mlon,mdepth,north,south,east,west,cont_int)
% This routine called to create base map used by rapid_widgit_lite

% Start a new figure if for web then no Matlab plot
if web_dis == 1
	figure('Visible','off')
else
	figure
end

% Navigation plot onhte left of the page will put text on the right later
subplot(1,2,1)

% Set up m_map plot
m_proj('lambert','lon',[west, east],'lat',[south, north]);
m_grid('box','on','color','k','linewidth',[1],'fontsize',[12]);
hold on

% Plot location of mooring
ha = m_plot(mlon,mlat,'k+');
set(ha,'MarkerSize',12,'LineWidth',1.5)

if strmatch(this_moor,'ship') == 1
    return
elseif web_dis == 0
	return
else
  saveas(gcf,[work_dir this_moor '_map'],'fig');
end
