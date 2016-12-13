close(figure(1));
close(figure(2));
fig1 = figure(1);
set(fig1, 'Position', [55 408 500 635],...
  'Color', [0.702 0.702 0.702], 'Numbertitle', 'off');

fig2 = figure(2);
set(fig2, 'Position', [575 448 600 540], 'Numbertitle', 'off');

regelkleur=[0.25 0.25 0.25];
knopkleur=[0.775 0.775 0.775];
editkleur=[0.9 0.9 0.9];
foutkleur=[0 0 0];

%******************************************************************

uicontrol(fig1,'Style','frame',...
          'Position',[0 625 510 5],'Background',regelkleur);

uicontrol(fig1, 'Style', 'text', 'Fontsize', 14,...
  'Position', [90 565 200 30],...
  'String', 'Diffusion constant k ',...
  'Horizontal','right');

setdif = uicontrol(fig1, 'Style', 'edit', 'Fontsize', 12,...
  'Position', [295 575 70 21],...
  'Background',knopkleur, ...
  'String', '0.01',...
  'Horizontal','left',...
  'Callback', ...
   ['figure(2),'...
    'kstr = get(setdif, ''String'');']);

kstr='0.01';

uicontrol(fig1,'Style','frame',...
          'Position',[0 545 510 5],'Background',regelkleur);

setroostertext = uicontrol(fig1, 'Style', 'text', 'Fontsize', 14,...
  'Position', [0 510 500 20],...
  'String', 'Grid');

setuniform = uicontrol(fig1, 'Style', 'checkbox', 'Fontsize', 12,...
  'Position', [100 470 100 30],...
  'Background',knopkleur, ...
  'String', 'uniform','Horizontal','center','Value',1,...
  'Callback',...
   [   'figure(2),'...
       'if (get(setupwindu,''Value'')==1),'...
          'methode=''Uuni''; '...
       'end,'...
       'if (get(setcentraal,''Value'')==1),'...
          'methode=''Cuni''; '...
       'end,'...
       'set(setuniform,''Value'',1),'...
       'set(setgerekt,''Value'',0),'...
       'set(setupwindu,''Visible'',''on''),'...
       'set(setcentraal,''Visible'',''on''),'...
       'set(setupwindg,''Visible'',''off''),'...
       'set(setcentraala,''Visible'',''off''),'...
       'set(setcentraalb,''Visible'',''off''),'...
       'rooster=''uni''; ' ]);

rooster='uni';

setgerekt = uicontrol(fig1, 'Style', 'checkbox', 'Fontsize', 12,...
  'Position', [300 470 100 30],...
  'Background',knopkleur, ...
  'String', 'stretched', 'Value',0,...
  'Callback',...
   [   'figure(2),'...
       'if (get(setupwindg,''Value'')==1),'...
          'methode=''Urek''; '...
       'end,'...
       'if (get(setcentraala,''Value'')==1),'...
          'methode=''Arek''; '...
       'end,'...
       'if (get(setcentraalb,''Value'')==1),'...
          'methode=''Brek''; '...
       'end,'...
       'set(setgerekt,''Value'',1),'...
       'set(setuniform,''Value'',0),'...
       'set(setupwindu,''Visible'',''off''),'...
       'set(setcentraal,''Visible'',''off''),'...
       'set(setupwindg,''Visible'',''on''),'...
       'set(setcentraala,''Visible'',''on''),'...
       'set(setcentraalb,''Visible'',''on''),'...
       'rooster=''rek''; ']);

%******************************************************************

uicontrol(fig1,'Style','frame',...
          'Position',[0 450 510 5],'Background',regelkleur);

setroostertext = uicontrol(fig1, 'Style', 'text', 'Fontsize', 14,...
  'Position', [10 410 475 20],...
  'String', 'Discretization');

%als uniform

setupwindu = uicontrol(fig1, 'Style', 'checkbox', 'Fontsize', 12,...
  'Position', [100 370 100 30],...
  'Background',knopkleur, ...
  'String', 'upwind', 'Visible','on','Value',1,...
  'Callback',...
   [   'figure(2),'...
       'set(setupwindu,''Value'',1),'...
       'set(setcentraal,''Value'',0),'...
       'methode=''Uuni''; ']);

setcentraal = uicontrol(fig1, 'Style', 'checkbox', 'Fontsize', 12,...
  'Position', [300 370 100 30],...
  'Background',knopkleur, ...
  'String', 'central', 'Visible','on','Value',0,...
  'Callback',...
   [   'figure(2),'...
       'set(setcentraal,''Value'',1),'...
       'set(setupwindu,''Value'',0),'...
       'methode=''Cuni''; ']);

methode='Uuni';

%als gerekt

setupwindg = uicontrol(fig1, 'Style', 'checkbox', 'Fontsize', 12,...
  'Position', [25 370 100 30],...
  'Background',knopkleur, ...
  'String', 'upwind','Visible','off','Value',1,...
  'Callback',...
   [   'figure(2),'...
       'set(setupwindg,''Value'',1),'...
       'set(setcentraala,''Value'',0),'...
       'set(setcentraalb,''Value'',0),'...
       'methode=''Urek''; ']);

setcentraala = uicontrol(fig1, 'Style', 'checkbox', 'Fontsize', 12,...
  'Position', [200 370 100 30],...
  'Background',knopkleur, ...
  'String', 'central A', 'Visible','off','Value',0,...
  'Callback',...
   [   'figure(2),'...
       'set(setupwindg,''Value'',0),'...
       'set(setcentraala,''Value'',1),'...
       'set(setcentraalb,''Value'',0),'...
       'methode=''Arek''; ']);

setcentraalb = uicontrol(fig1, 'Style', 'checkbox', 'Fontsize', 12,...
  'Position', [375 370 100 30],...
  'Background',knopkleur, ...
  'String', 'central B', 'Visible','off','Value',0,...
  'Callback',...
   [   'figure(2),'...
       'set(setupwindg,''Value'',0),'...
       'set(setcentraala,''Value'',0),'...
       'set(setcentraalb,''Value'',1),'...
       'methode=''Brek''; ']);

uicontrol(fig1,'Style','frame',...
          'Position',[0 350 510 5],'Background',regelkleur);

%******************************************************************

setitertext = uicontrol(fig1, 'Style', 'text', 'Fontsize', 14,...
  'Position', [10 315 475 20],...
  'String', 'Iterative method');

setjor = uicontrol(fig1, 'Style', 'checkbox', 'Fontsize', 12,...
  'Position', [100 285 100 30],...
  'Background',knopkleur, ...
  'String', 'JOR','Value',1,...
  'Callback', ...
    [ 'figure(2),'...
      'set(setjor,''Value'',1),'...
      'set(setsor,''Value'',0),'...
      'itmethode=''JOR''; ']);

setsor = uicontrol(fig1, 'Style', 'checkbox', 'Fontsize', 12,...
  'Position', [300 285 100 30],...
  'String', 'SOR','Value',0,...
  'Background',knopkleur, ...
  'Callback', ...
    [ 'figure(2),'...
      'set(setjor,''Value'',0),'...
      'set(setsor,''Value'',1),'...
      'itmethode=''SOR''; ']);

itmethode='JOR';

uicontrol(fig1, 'Style', 'text', 'Fontsize', 12,...
  'Position', [180 235 120 30],...
  'String', 'relaxation factor w ',...
  'Horizontal','right');

setpar = uicontrol(fig1, 'Style', 'edit', 'Fontsize', 12,...
  'Position', [220 214 50 21],...
  'Background',knopkleur, ...
  'String', '1',...
  'Horizontal','left',...
  'Callback', ...
   ['figure(2),'...
    'wstr = get(setpar, ''String'');']);

wstr='1';

%******************************************************************

uicontrol(fig1,'Style','frame',...
          'Position',[0 190 510 5],'Background',regelkleur);

pushprog = uicontrol(fig1, 'Style', 'pushbutton', 'Fontsize', 16,...
  'Position', [0 144 510 35],...
  'Background',[0.0 0.50 0.0], ...
  'Foreground',[1 1 1], ...
  'String', 'EXECUTE PROGRAM',...
  'Callback', 'figure(2),ns_iterative');

pushexit = uicontrol(fig1, 'Style', 'pushbutton', 'Fontsize', 16,...
  'Position', [0 95 510 35],...
  'Background',[0.50 0.0 0.0], ...
  'Foreground',[1 1 1], ...
  'String', 'STOP PROGRAM EXECUTION',...
  'Callback','figure(2),keyboard');

uicontrol(fig1,'Style','frame',...
          'Position',[0 76 510 5],'Background',regelkleur);

%******************************************************************

pushexit = uicontrol(fig1, 'Style', 'pushbutton', 'Fontsize', 16,...
  'Position', [0 27 510 35],...
  'Background',foutkleur, ...
  'Foreground',[1 1 1], ...
  'String', 'CLOSE CONTROL PANEL',...
  'Callback', ['close(fig2); close(fig1)']);

uicontrol(fig1,'Style','frame',...
          'Position',[0 8 510 5],'Background',regelkleur);