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

C = cellfun(@(x,y,z) [x y z], resultsallU8_1_, resultsallU8_2_, resultsallU8_3_, 'UniformOutput', false);
C = cellfun(@(x,y,z) [x y z], resultsallU8_1, resultsallU8_4, resultsallU8_3, 'UniformOutput', false);

Gall = [0.001,0.01,0.1,0.2,0.4,0.5,0.6,0.7,0.9,1.0,1.2,1.7,2,2.5,3,4,5,6,7,8,9,10,11,15,20,25];


disp(C{1})   % → 1×8 cell: {'a1','a2','a3','a4','f1','f2','f3','f4'}
