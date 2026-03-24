A = { {'a1','a2','a3','a4'}; 
      {'b1','b2','b3','b4'}; 
      {'c1','c2','c3','c4'}; 
      {'d1','d2','d3','d4'}; 
      {'e1','e2','e3','e4'} };

B = { {'f1','f2','f3','f4'}; 
      {'g1','g2','g3','g4'}; 
      {'h1','h2','h3','h4'}; 
      {'i1','i2','i3','i4'}; 
      {'j1','j2','j3','j4'} };

C = cellfun(@(x,y,z,w) [x y z w], ...
            resultsallU8_4_2, resultsallU8_f6, resultsallU8_f11, resultsallF17, ...
            'UniformOutput', false);

C = cellfun(@(x,y,z) [x y z], resultsallU8_1, resultsallU8_4, resultsallU8_3, 'UniformOutput', false);

Gall = [0.001,0.01,0.1,0.2,0.4,0.5,0.6,0.7,0.9,1.0,1.2,1.7,2,2.5,3,4,5,6,7,8,9,10,11,15,20,25];


disp(C{1})   % → 1×8 cell: {'a1','a2','a3','a4','f1','f2','f3','f4'}


             % 14×1 cell, each cell is 1×2
C1 = reshape(resultsf23.', 1, []); 

C1 = num2cell([resultsf23{:}]);

resultsC = [];
for i = 1:28
    resutlsC = [resultsC, C1{i}];
end

% C14x1: 14×1 cell, each entry is (supposedly) a 1×2 container

tmp = cellfun(@(x) x(:).', resultsMar24, 'UniformOutput', false);  % each -> 1×2 (still could be cell/struct)
tmp = [tmp{:}];                                             % -> 1×28 (cell or numeric)
C4 = {tmp};

results = [C4];


