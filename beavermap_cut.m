%% Configuration
ax = [729886.11,730486.35,4701464.75, 4701796.03]; % domain
step_vine = 2.25; % distance b/w vines
coord = readtable("HARV_Model_Trees.csv");
tree_x = coord.average_easting;
tree_y = coord.average_northing;
tree_h = coord.height;
tree_n = length(tree_x);

xPoly = [730246.10, 730213.50, 730143.03, 730069.58, 730043.36, 730083.36, 730111.31, 730177.93];
yPoly = [4701622.69, 4701664.86, 4701689.10, 4701701.14, 4701689.71, 4701603.09, 4701576.89, 4701559.74];

indexescut = 1:tree_n ;
donotcut = @(i) not(inpolygon(tree_x(i), tree_y(i), xPoly, yPoly));
notcut = donotcut(indexescut)
tree_x = tree_x(notcut);
tree_y = tree_y(notcut);
tree_h = tree_h(notcut);
tree_n = length(tree_x);

%Adding the vineyard patch
[tree_x2, tree_y2] = vinemaker(730126.50,	4701673.33,	730149.78,	4701678.99,	730145.04,	4701605.57,	730171.43,	4701614.45);
tree_n2 = length(tree_x2);
tree_h2 = zeros(1,tree_n2) + 2.5;

tree_x = [tree_x ; tree_x2'];
tree_y = [tree_y ; tree_y2'];
tree_h = [tree_h ; tree_h2'];
tree_n = tree_n + tree_n2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vinexs, vineys] = vinemaker(TopLeftE,	TopLeftN,	TopRightE,	TopRightN,	BottLeftE,	BottLeftN,	BottRightE,	BottRightN)
    h = 2.25; % Vine step size;
    vinexs = [];
    vineys = [];

    theta = atan((TopRightN-TopLeftN)/(TopRightE-TopLeftE));
    alpha = atan((BottLeftE-TopLeftE)/(BottLeftN-TopLeftN));
    nrow = floor((TopRightE-TopLeftE)/(h * cos(theta))) + 1;
    row1xs = linspace(TopLeftE, TopRightE, nrow);
    row1ys = linspace(TopLeftN, TopRightN, nrow);
    while ~(any((row1xs > BottRightE) & (row1ys < BottRightN)) || any((row1xs < BottLeftE) & (row1ys < BottLeftN)))
        vinexs = [vinexs row1xs];
        vineys = [vineys row1ys];
        row1xs = row1xs - h*sin(alpha);
        row1ys = row1ys - h*cos(alpha);
    end

end