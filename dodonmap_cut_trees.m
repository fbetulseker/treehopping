%% Configuration
ax = [359084.89,359990.42,4306816.41, 4307379.35]; % domain
step_vine = 2.25; % distance b/w vines
coord = readtable("SERC_Model_Trees.csv");
tree_x = coord.average_easting;
tree_y = coord.average_northing;
tree_h = coord.height;
tree_n = length(tree_x);
nvines = 10000;

xPoly = [359240.41, 359425.47, 359716.36, 359778.23, 359337.78];
yPoly = [4307134.57, 4307331.54, 4307348.36, 4307168.75, 4306947.00];

indexescut = 1:tree_n ;
donotcut = @(i) not(inpolygon(tree_x(i), tree_y(i), xPoly, yPoly));
notcut = donotcut(indexescut);
tree_x = tree_x(notcut);
tree_y = tree_y(notcut);
tree_h = tree_h(notcut);
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

tree_x = tree_x-min(tree_x);
tree_y =tree_y-min(tree_y);

tree_h = [tree_h ; tree_h2'];
tree_n = tree_n + tree_n2;



% Adding a vineyard patch 6 - triangle under patch 3
[tree_x2, tree_y2] = vinemaker_triangle3(359557.98,	4307134.14,	359764,	4307170.65, 359685.5, 4307119.78, nvines);
tree_n2 = length(tree_x2);
tree_h2 = zeros(1,tree_n2) + 2.5;

tree_x = [tree_x ; tree_x2'-ax(1)];
tree_y = [tree_y ; tree_y2'-ax(3)];
tree_h = [tree_h ; tree_h2'];
tree_n = tree_n + tree_n2;

ax = ax-[ax(1),ax(1),ax(3),ax(3)]; % domain

border_trees = tree_y > (ax(4)-30);



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


function [vinexs, vineys] = vinemaker_empty(TopLeftE,	TopLeftN,	TopRightE,	TopRightN,	BottLeftE,	BottLeftN,	BottRightE,	BottRightN, nvines)
    h = 2.25; % Vine step size;
    vinexs = [];
    vineys = [];

    theta = atan((TopRightN-TopLeftN)/(TopRightE-TopLeftE));
    alpha = atan((BottLeftE-TopLeftE)/(BottLeftN-TopLeftN));
    nrow = floor((TopRightE-TopLeftE)/(h * cos(theta))) + 1;
    row1xs = linspace(TopLeftE, TopRightE, nrow);
    row1ys = linspace(TopLeftN, TopRightN, nrow);
    if length(row1xs) > 2*nvines
        row1xs_empty = row1xs([1:nvines, length(row1xs)+1-nvines:length(row1xs)]);
        row1ys_empty = row1ys([1:nvines, length(row1ys)+1-nvines:length(row1ys)]);
    else
        row1xs_empty = row1xs;
        row1ys_empty = row1ys;
    end
    
    rowcounter = 1;
    while ~(any((row1xs > BottRightE) & (row1ys < BottRightN)) || any((row1xs < BottLeftE) & (row1ys < BottLeftN)))
        nfuturerowxs = row1xs - (nvines)*h*sin(alpha);
        nfuturerowys = row1ys - (nvines)*h*cos(alpha);
        if ~(any((nfuturerowxs > BottRightE) & (nfuturerowys < BottRightN)) || any((nfuturerowxs < BottLeftE) & (nfuturerowys < BottLeftN))) && rowcounter > nvines
            vinexs = [vinexs row1xs_empty];
            vineys = [vineys row1ys_empty];
        else
            vinexs = [vinexs row1xs];
            vineys = [vineys row1ys];
        end

        row1xs = row1xs - h*sin(alpha);
        row1ys = row1ys - h*cos(alpha);
        row1xs_empty = row1xs_empty - h*sin(alpha);
        row1ys_empty = row1ys_empty - h*cos(alpha);
        rowcounter = rowcounter + 1;

    end

end

function [vinexs, vineys] = vinemaker_empty3(TopLeftE,	TopLeftN,	TopRightE,	TopRightN,	BottLeftE,	BottLeftN,	BottRightE,	BottRightN, nvines)
    h = 2.25; % Vine step size;
    vinexs = [];
    vineys = [];

    theta = atan((TopRightN-TopLeftN)/(TopRightE-TopLeftE));
    alpha = atan((BottLeftE-TopLeftE)/(BottLeftN-TopLeftN));
    nrow = floor((TopRightE-TopLeftE)/(h * cos(theta))) + 1;
    row1xs = linspace(TopLeftE, TopRightE, nrow);
    row1ys = linspace(TopLeftN, TopRightN, nrow);
    if length(row1xs) > 2*nvines
        row1xs_empty = row1xs([1:nvines, length(row1xs)+1-nvines:length(row1xs)]);
        row1ys_empty = row1ys([1:nvines, length(row1ys)+1-nvines:length(row1ys)]);
    else
        row1xs_empty = row1xs;
        row1ys_empty = row1ys;
    end
    
    rowcounter = 1;
    while ~(any((row1xs > BottRightE) & (row1ys < BottRightN)) || any((row1xs < BottLeftE) & (row1ys < BottLeftN)))
        nfuturerowxs = row1xs - (nvines)*h*sin(alpha);
        nfuturerowys = row1ys - (nvines)*h*cos(alpha);
        if rowcounter > nvines
            vinexs = [vinexs row1xs_empty];
            vineys = [vineys row1ys_empty];
        else
            vinexs = [vinexs row1xs];
            vineys = [vineys row1ys];
        end

        row1xs = row1xs - h*sin(alpha);
        row1ys = row1ys - h*cos(alpha);
        row1xs_empty = row1xs_empty - h*sin(alpha);
        row1ys_empty = row1ys_empty - h*cos(alpha);
        rowcounter = rowcounter + 1;

    end

end

function [vinexs, vineys] = vinemaker_triangle(TopLeftE,	TopLeftN,	TopRightE,	TopRightN,	BottLeftE,	BottLeftN, nvines)
    h = 2.25; % Vine step size;
    vinexs = [];
    vineys = [];

    theta = atan((TopRightN-TopLeftN)/(TopRightE-TopLeftE));
    alpha = atan((BottLeftE-TopLeftE)/(BottLeftN-TopLeftN));
    rowcounter = 1;
    while TopLeftE < BottLeftE && TopLeftN > BottLeftN
        tmpx = TopLeftE;
        tmpy = TopLeftN;
        tmprowx = [];
        tmprowy = [];
        while tmpx < (tmpy-BottLeftN)*(TopRightE-BottLeftE)/(TopRightN-BottLeftN)+BottLeftE
            tmprowx = [tmprowx tmpx];
            tmprowy = [tmprowy tmpy];
            tmpx = tmpx + h*cos(theta);
            tmpy = tmpy + h*sin(theta);
        end
        
        if rowcounter > nvines && length(tmprowx) > 2*nvines
            tmprowx = tmprowx([1:nvines, length(tmprowx)+1-nvines:length(tmprowx)]);
            tmprowy = tmprowy([1:nvines, length(tmprowy)+1-nvines:length(tmprowy)]);
        end

        vinexs = [vinexs tmprowx];
        vineys = [vineys tmprowy];

        TopLeftE = TopLeftE + h*sin(alpha)/cos(theta-alpha);
        TopLeftN = TopLeftN - h*cos(alpha)/cos(theta-alpha);
        rowcounter = rowcounter + 1;
    end

end

function [vinexs, vineys] = vinemaker_triangle3(TopLeftE,	TopLeftN,	TopRightE,	TopRightN,	BottLeftE,	BottLeftN, nvines)
    h = 2.25; % Vine step size;
    vinexs = [];
    vineys = [];
    

    theta = atan((TopRightN-TopLeftN)/(TopRightE-TopLeftE));
    alpha = atan((BottLeftE-TopLeftE)/(TopLeftN-BottLeftN));

    TopLeftE = TopLeftE ;
    TopLeftN = TopLeftN;
    rowcounter = 1;
    while TopLeftE < BottLeftE && TopLeftN > BottLeftN
        tmpx = TopLeftE;
        tmpy = TopLeftN;
        tmprowx = [];
        tmprowy = [];
        while tmpx < (tmpy-BottLeftN)*(TopRightE-BottLeftE)/(TopRightN-BottLeftN)+BottLeftE
            tmprowx = [tmprowx tmpx];
            tmprowy = [tmprowy tmpy];
            tmpx = tmpx + h*cos(theta);
            tmpy = tmpy + h*sin(theta);
        end
        
        if length(tmprowx) > 2*nvines
            tmprowx = tmprowx([1:nvines, length(tmprowx)+1-nvines:length(tmprowx)]);
            tmprowy = tmprowy([1:nvines, length(tmprowy)+1-nvines:length(tmprowy)]);
        end

        vinexs = [vinexs tmprowx];
        vineys = [vineys tmprowy];

        TopLeftE = TopLeftE +  h*sin(alpha)/cos(theta-alpha);
        TopLeftN = TopLeftN - h*cos(alpha)/cos(theta-alpha);
        rowcounter = rowcounter + 1;
    end

end