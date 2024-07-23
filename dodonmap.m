%% Configuration
ax = [359084.89,359990.42,4306816.41, 4307379.35]; % domain
step_vine = 2.25; % distance b/w vines
coord = readtable("SERC_Model_Trees.csv");
tree_x = coord.average_easting;
tree_y = coord.average_northing;
tree_h = coord.height;
tree_n = length(tree_x);


% Adding a vineyard patch 4
[tree_x2, tree_y2] = vinemaker(359464.57,	4307115.71	,359508.15	,4307132.25	,359527.11	,4306988.37	,359572.47	,4307017.21);
tree_n2 = length(tree_x2);
tree_h2 = zeros(1,tree_n2) + 2.5;

tree_x = [tree_x ; tree_x2'];
tree_y = [tree_y ; tree_y2'];
tree_h = [tree_h ; tree_h2'];
tree_n = tree_n + tree_n2;


% Adding a vineyard patch 1
[tree_x2, tree_y2] = vinemaker(359267.76,	4307160.5,	359396.96,	4307216.37,	359293.06,	4307105.1,	359423.29,	4307161.06);
tree_n2 = length(tree_x2);
tree_h2 = zeros(1,tree_n2) + 2.5;

tree_x = [tree_x ; tree_x2'];
tree_y = [tree_y ; tree_y2'];
tree_h = [tree_h ; tree_h2'];
tree_n = tree_n + tree_n2;


% Adding a vineyard patch 2
[tree_x2, tree_y2] = vinemaker(359527.87,	4307329.39,	359697.66,	4307340.81,	359539.92,	4307231.04,	359703.17,	4307260.23);
tree_n2 = length(tree_x2);
tree_h2 = zeros(1,tree_n2) + 2.5;

tree_x = [tree_x ; tree_x2'];
tree_y = [tree_y ; tree_y2'];
tree_h = [tree_h ; tree_h2'];
tree_n = tree_n + tree_n2;


% Adding a vineyard patch 3
[tree_x2, tree_y2] = vinemaker(359544.66,	4307214.86,	359750.34,	4307254.17,	359557.98,	4307134.14,	359759.09,	4307170.65);

xmax = 359759.09;
xmin = 359557.98;
ymax = 4307170.65;
ymin = 4307119.78;

treecheck = find(tree_x < xmax + 20 & tree_x > xmin - 20 & tree_y < ymax + 20 & tree_y > ymin - 40);
treecheckx = tree_x(treecheck);
treechecky = tree_y(treecheck);
numbertrial = 0;


tree_n2 = length(tree_x2);
tree_h2 = zeros(1,tree_n2) + 2.5;

tree_x = [tree_x ; tree_x2'];
tree_y = [tree_y ; tree_y2'];
tree_h = [tree_h ; tree_h2'];
tree_n = tree_n + tree_n2;



% Adding a vineyard patch 5
[tree_x2, tree_y2] = vinemaker(359533.39,	4306975.83,	359602.82,	4307016.89,	359582.08,	4306869.83,	359659.14,	4306920.31);
tree_n2 = length(tree_x2);
tree_h2 = zeros(1,tree_n2) + 2.5;

tree_x = [tree_x ; tree_x2'];
tree_y = [tree_y ; tree_y2'];
tree_h = [tree_h ; tree_h2'];
tree_n = tree_n + tree_n2;



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