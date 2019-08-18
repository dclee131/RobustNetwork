function VisualizeNetwork(sys_case,NodalValues,FlowValues,node_scale,flow_scale,ref_bus)
height = 0;
hold all; box on;
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);

if sys_case==9
    NodeCoordinates=[40 5; 0 40; 80 40; 40 10; 5 40; 75 40; 40 20; 30 30; 50 30; 20 40; 40 40; 60 40];
    grid_edge=[9 8; 7 8; 9 6; 7 5; 5 4; 6 4; 2 7; 3 9; 1 4];
    idx_gen=[1 2 3]'; num_gen=size(idx_gen,1);
    EdgeConnections=[[1:num_gen]' num_gen+idx_gen; num_gen+grid_edge];
    axis([-5 85 0 50])
elseif sys_case==14
    gen_node=[10 55; 75 8; 15 5; 90 75; 25 50];
    NodeCoordinates=[gen_node; 10 50; 15 10; 75 13; 75 55; 20 35; 25 55; 82 66; 90 70; 72 72; 50 70; 38 75; 15 85; 30 90; 50 85];
    grid_edge=[2 5; 6 12;12 13; 6 13; 6 11; 11 10; 9 10; 9 14; 14 13; 7 9; 1 2; 3 2; 3 4; 1 5; 5 4; 2 4; 5 6; 4 9; 4 7; 8 7];
    idx_gen=[1 3 2 8 6]'; num_gen=size(idx_gen,1);
    EdgeConnections=[[1:num_gen]' num_gen+idx_gen; num_gen+grid_edge];
elseif sys_case==39
    gen_node=[20 65; 30 5; 45 5; 70 5; 60 5; 90 35; 80 5; 40 105; 80 70; 30 105];
    NodeCoordinates=[gen_node; 25 70; 30 80; 30 70; 35 60; 30 50; 30 20; 25 30; 25 40; 20 50; 45 25; 35 35; 45 45; 55 35; 55 60;
    62.5 65; 70 70; 50 70; 40 70; 70 40; 60 25; 85 60; 85 50; 80 20; 75 50; 40 90; 50 90; 50 80; 65 80; 80 90; 30 100;
    30 10; 45 10; 70 10; 60 10; 90 40; 80 10; 40 100; 80 75; 20 60];
    grid_edge=[1 2; 1 39; 2 3; 2 25; 2 30; 3 4; 3 18; 4 5; 4 14; 5 8; 6 5; 6 7; 6 11; 7 8; 8 9; 9 39; 
    10 11; 10 13; 10 32; 12 11; 12 13; 13 14; 14 15; 15 16; 16 17; 16 19; 16 21; 16 24; 17 18; 17 27; 19 33;
    19 20; 20 34; 21 22; 22 23; 22 35; 23 24; 23 36; 25 26; 25 37; 26 27; 26 28; 26 29; 28 29; 29 38; 6 31];
    idx_gen=[39 31 32 33 34 35 36 37 38 30]'; num_gen=size(idx_gen,1);
    EdgeConnections=[[1:num_gen]' num_gen+idx_gen; num_gen+grid_edge];
    axis([15 95 0 110])
end

for i = 1:size(FlowValues,1)
   nodeA = [NodeCoordinates(EdgeConnections(i,1),1) NodeCoordinates(EdgeConnections(i,1),2)];
   nodeB = [NodeCoordinates(EdgeConnections(i,2),1) NodeCoordinates(EdgeConnections(i,2),2)];
   thickness = abs(FlowValues(i))*flow_scale(2)+flow_scale(1); % [0 10]
   line([nodeA(1) nodeB(1)],[nodeA(2) nodeB(2)],[height   height], 'LineWidth', thickness, 'color', 'k','LineStyle','-');
end

t=0:2*pi/20:2*pi;
for i = 1:size(NodalValues,1)
   center = [NodeCoordinates(i,1) NodeCoordinates(i,2)];
   radius = abs(NodalValues(i))*node_scale(2)+node_scale(1); % [0.2 0.5]
   if i==ref_bus
       fill(sin(t)*radius+center(1),cos(t)*radius+center(2),'r','LineWidth',1);
   elseif i<=num_gen
       fill(sin(t)*radius+center(1),cos(t)*radius+center(2),'b','LineWidth',1);
   else
       fill(sin(t)*radius+center(1),cos(t)*radius+center(2),'w','LineWidth',1);
   end
end


