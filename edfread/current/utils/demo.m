% sample usage of utility functions

% read data w/ edfread
[dat, meta] = edfread('/Users/jsteger/lagTest/AB_VP016.EDF')

% lets do some pre-processing:
% convert the experimentor's name to a numerical id
experimentors = {'johannes'; 'andreas'; 'dummy'}
[discard, discard, index] = intersect(meta.EXPERIMENTOR, experimentors);
% (index is now the index of the string 'meta.EXPERIMENTOR' in the cell
% array 'experimentors', e.g. 2 for 'andreas')
% note: fixations.m expects the string representation of an integer here,
% so we convert the index
meta.EXPERIMENTORS = num2str(index);

% now we convert to fixmat
fixmat0 = fixations(dat, meta);

% say we did this for two other subjects too and saved the results to
% fixmat1 and fixmat2 (in practice, you will of course want to do that
% within a function) 
% we can now merge all data into one structure:
fixmat = concat(fixmat0, fixmat1, fixmat2);

% you should really save that fixmat now
save('precious.mat', 'fixmat');

% some use cases
% plot all fixations condition 0 and 1 with different colors
condfix = split(fixmat, 'condition', [0, 1]);

plot(condfix(1).x, condfix(1).y, '+r');
hold;
plot(condfix(2).x, condfix(2).y, '+b');