% This code read the .mat files generate from python
% and then calculate null space of Sint and generate input files for GAMS

%if file is in another directory
filepath1 = './optstoic_db_v3/data_1/';
INPUTFILE = fullfile(filepath1, 'Sint_no_cofactor_20160831.mat');
load(INPUTFILE);
% OUTPUTF1 = fullfile(filepath, '20160622_gams_null_sij_nocofactor.txt');
% OUTPUTF2 = fullfile(filepath, '20160622_gams_loops_nocofactor.txt');
filepath = './optstoic_db_v3/data_1/';
fprefix = 'optstoic_v3';
OUTPUTF1 = fullfile(filepath, strcat(fprefix, '_null_sij_nocofactor_20170627.txt'));
OUTPUTF2 = fullfile(filepath, strcat(fprefix, '_loops_nocofactor_20170627.txt'));

if (exist('Sint_sparse') ~= 1)
    disp('Variable Sint_sparse is not available.')
end

if (exist('reactionList') ~= 1)
    disp('Variable reactionList is not available.')    
end


fprintf('Calculating the null...\n');

Sint = spconvert(Sint_sparse);
%Get the rational basis of the null space of Sint
Nint = null(Sint,'r');

%%1. Clean up the Nint with values < 1e-9
eps = 1e-9;
indices = find(abs(Nint)<eps);
Nint(indices) = 0;

%Remove single reaction loop (reaction involving only cofactors)
singleReactionLoop = [];
for col=1:size(Nint, 2)
    countRxn = find(abs(Nint(:,col))>1e-9);
    if length(countRxn) == 1
        fprintf('single reaction loop found! %s\n',reactionList{col})
        singleReactionLoop = [singleReactionLoop; col];
        Nint(countRxn, col) = 0;
    end
end
singleReaction = reactionList(singleReactionLoop);

%%writing input for GAMS
%%the Nint matrix
fprintf('Writing output files...\n');
%%2. Get indices and values for all nonzero elements
[jVector, kVector, val] = find(Nint);
%%Herein, j is reaction; k is pathway/loop

%%3. Write to GAMS include file
fileID2=fopen(OUTPUTF1, 'w');
for n = 1:length(jVector)
    fprintf(fileID2,'''L%d''.''%s'' %f\n',kVector(n), reactionList{jVector(n)}, val(n));
end
fclose(fileID2);

%%4. Write the elements for the set loops
fileID3=fopen(OUTPUTF2, 'w');
fprintf(fileID3,'''L%d''\n', 1:size(Nint,2));
fclose(fileID3);