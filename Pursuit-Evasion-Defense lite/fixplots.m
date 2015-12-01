function h =  fixplots(f)
% This function works in MATLAB 2014b and later versions.
% Open a figure and manually hf to it by doing hf = gcf.

h = f.Children;
f.Position = f.Position+[0 0 0 100];

for i=1:4
    h(i).Units = 'pixels';
    h(i).ActivePositionProperty = 'outerposition';
    h(i).Title.String='';
    h(i).XTick=[-1  0  1];
    h(i).YTick=[-1  0  1];
    if i>=3,h(i).XTick = [0 1];end;
    h(i).Box = 'on';
    h(i).Position = h(i).Position + [-15 85 30 30] + [0 -10 0 0]*(i<=2);
    
    hp = h(i).Children; 
    if length(hp)==7
        %Default order: Basin,Obstacle,Evader,Capture,Interception,Defender,Reach
        %New order:     Evader,Defender,Capture,Interception,Obstacle,Reach,Basin
    end
    % Figures for full game are fine by default
    h(i).FontName='Helvetica';
    h(i).FontSize=12;
end