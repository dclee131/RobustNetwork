function VisualizeNetwork(sys_case,NodalValues,FlowValues)
height = 0;
hold on; box on;
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);

if sys_case==9
    NodeCoordinates=[40 10; 5 40; 75 40; 40 20; 30 30; 50 30; 20 40; 40 40; 60 40];
    EdgeConnections=[ 9 8; 7 8; 9 6; 7 5; 5 4; 6 4; 2 7; 3 9; 1 4];
    axis([0 80 0 50])
elseif sys_case==14
    NodeCoordinates=[10 50; 15 10; 75 13; 75 55; 20 35; 25 55; 82 66; 90 70; 72 72; 50 70; 38 75; 15 85; 30 90; 50 85];
    EdgeConnections=[2 5; 6 12;12 13; 6 13; 6 11; 11 10; 9 10; 9 14; 14 13; 7 9; 1 2; 3 2; 3 4; 1 5; 5 4; 2 4; 5 6; 4 9; 4 7; 8 7];
elseif sys_case==39
    NodeCoordinates=[25 70; 30 80; 30 70; 35 60; 30 50; 30 20; 25 30; 25 40; 20 50; 45 25; 35 35; 45 45; 55 35; 55 60;
    62.5 65; 70 70; 50 70; 40 70; 70 40; 60 25; 85 60; 85 50; 80 20; 75 50; 40 90; 50 90; 50 80; 65 80; 80 90; 30 100;
    30 10; 45 10; 70 10; 60 10; 90 40; 80 10; 40 100; 80 75; 20 60];
    EdgeConnections=[1 2; 1 39; 2 3; 2 25; 2 30; 3 4; 3 18; 4 5; 4 14; 5 8; 6 5; 6 7; 6 11; 7 8; 8 9; 9 39; 
    10 11; 10 13; 10 32; 12 11; 12 13; 13 14; 14 15; 15 16; 16 17; 16 19; 16 21; 16 24; 17 18; 17 27; 19 33;
    19 20; 20 34; 21 22; 22 23; 22 35; 23 24; 23 36; 25 26; 25 37; 26 27; 26 28; 26 29; 28 29; 29 38; 6 31];
    axis([10 95 5 105])
end

for i = 1:size(FlowValues,1)
   nodeA = [NodeCoordinates(EdgeConnections(i,1),1) NodeCoordinates(EdgeConnections(i,1),2)];
   nodeB = [NodeCoordinates(EdgeConnections(i,2),1) NodeCoordinates(EdgeConnections(i,2),2)];
   thickness = abs(FlowValues(i))*10;
   if FlowValues(i)==0
       line([nodeA(1) nodeB(1)],[nodeA(2) nodeB(2)],[height   height], 'LineWidth', 1, 'color', 'r','LineStyle',':');
   elseif FlowValues(i)>0
       line([nodeA(1) nodeB(1)],[nodeA(2) nodeB(2)],[height   height], 'LineWidth', thickness, 'color', 'r','LineStyle','-');
       %quiver(nodeB(1), nodeB(2),(nodeA(1)-nodeB(1))/2, (nodeA(2)-nodeB(2))/2, 'LineWidth',  thickness, 'color', 'r','MaxHeadSize',300,'AutoScale','on');
   elseif FlowValues(i)<0
       line([nodeA(1) nodeB(1)],[nodeA(2) nodeB(2)],[height   height], 'LineWidth', thickness, 'color', 'r','LineStyle','-');
       %quiver(nodeA(1), nodeA(2),(nodeB(1)-nodeA(1))/2, (nodeB(2)-nodeA(2))/2, 'LineWidth', thickness, 'color', 'r','MaxHeadSize',300,'AutoScale','on');
   end
end

t=0:2*pi/20:2*pi;
for i = 1:size(NodalValues,1)
   center = [NodeCoordinates(i,1) NodeCoordinates(i,2)];
   radius = abs(NodalValues(i))*1+0.5;
   if NodalValues(i) > 0
       fill3(sin(t)*radius+center(1),cos(t)*radius+center(2),zeros(size(t))+height,'b');
   elseif NodalValues(i) < 0
       fill3(sin(t)*radius+center(1),cos(t)*radius+center(2),zeros(size(t))+height,'r');
   else
       fill3(sin(t)*radius+center(1),cos(t)*radius+center(2),zeros(size(t))+height,'w');
   end
end


