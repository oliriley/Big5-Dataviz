% Visualizing data from a set of >19k Big Five Inventory responses from an
% anonymized dataset:
% https://openpsychometrics.org/_rawdata/BIG5
close all
clear all

load BIG5_Data.mat
% Expected result: 19599x5 matrix of doubles.
% 
% Each row of the matrix is five numbers representing the sums of a single 
% respondant's answers to the questions of a 50-item Big 5 Inventory. Each 
% of five traits were assessed by 10 questions on a scale from 1-5, so
% values in the matrix will range from 10 (a respondant giving an answer of
% '1' to all questions relevant to a given trait) to 50 (giving an answer
% of '5' to all questions for a trait.

%%

figure
subplot(2,3,1)
plot(data(:,1))
subplot(2,3,2)
plot(data(:,2))
subplot(2,3,3)
plot(data(:,3))
subplot(2,3,4)
plot(data(:,4))
subplot(2,3,5)
plot(data(:,5))

% Examining the data shows one error response with 0 across the board.
% There are also a number of responses in which a respondant answered '1'
% or '5' for ALL questions, which strongly suggests they were not taking
% the questionnaire seriously. These have also been removed using:

loop = 1;
while loop <= length(data(:,1))
    if (sum(data(loop,:)) <= 50 || (sum(data(loop,:)) >= 250))
        data = [data(1:loop-1,:) ; data(loop+1:end,:)];
        loop = loop-1;
    end
    loop = loop +1;
end

% Some anomalous answers remain, but it's good enough for a first pass. 

% One optional modification: subtracting off the mean of each column and
% dividing the remaining columns by their variance.
% I believe, but have not yet checked, that this will make some of the later
% analysis easier, but haven't yet checked.)

% Subtract the means:
data = data - repmat([mean(data(:,1)),mean(data(:,2)),mean(data(:,3)),mean(data(:,4)),mean(data(:,5))],size(data,1),1);

% Divide by variance:
for loop = 1:size(data,2)
    data(:,loop) = data(:,loop)./std(data(:,loop));
end

figure
subplot(2,3,1)
plot(data(:,1))
subplot(2,3,2)
plot(data(:,2))
subplot(2,3,3)
plot(data(:,3))
subplot(2,3,4)
plot(data(:,4))
subplot(2,3,5)
plot(data(:,5))

%%

figure
subplot(2,3,1)
hist(data(:,1),[floor(min(data(:,1))):ceil(max(data(:,1)))])
title('Openness')

subplot(2,3,2)
hist(data(:,2),[floor(min(data(:,1))):ceil(max(data(:,2)))])
title('Conscientiousness')

subplot(2,3,3)
hist(data(:,3),[floor(min(data(:,1))):ceil(max(data(:,1)))])
title('Extraversion')

subplot(2,3,4)
hist(data(:,4),[floor(min(data(:,1))):ceil(max(data(:,1)))])
title('Agreeableness')

subplot(2,3,5)
hist(data(:,5),[floor(min(data(:,1))):ceil(max(data(:,1)))])
title('Neuroticism')

% Visualizing the histograms of each trait column individually produces
% Gaussian curves, suggesting the traits are normally distributed.

%%
figure
subplot(3,4,1)
plot(data(:,1),data(:,2),'.')
xlabel('O')
ylabel('C')

subplot(3,4,2)
plot(data(:,1),data(:,3),'.')
xlabel('O')
ylabel('E')

subplot(3,4,3)
plot(data(:,1),data(:,4),'.')
xlabel('O')
ylabel('A')

subplot(3,4,4)
plot(data(:,1),data(:,5),'.')
xlabel('O')
ylabel('N')

subplot(3,4,5)
plot(data(:,2),data(:,3),'.')
xlabel('C')
ylabel('E')

subplot(3,4,6)
plot(data(:,2),data(:,4),'.')
xlabel('C')
ylabel('A')

subplot(3,4,7)
plot(data(:,2),data(:,5),'.')
xlabel('C')
ylabel('N')

subplot(3,4,8)
plot(data(:,3),data(:,4),'.')
xlabel('E')
ylabel('A')

subplot(3,4,9)
plot(data(:,3),data(:,5),'.')
xlabel('E')
ylabel('N')

subplot(3,4,10)
plot(data(:,4),data(:,5),'.')
xlabel('A')
ylabel('N')

% Simply plotting each trait against the others is not very illustrative,
% because there is no many indication of how many datapoints lie at each
% intersection.

%%
% Instead, construct 3D histograms:
 
RANGE = [floor(min(min(data))):ceil(max(max(data)))];
figure
hist3([data(:,1),data(:,2)],{RANGE,RANGE})
title('Openness vs Concientiousness')
xlabel('O')
ylabel('C')
set(gcf,'renderer','opengl')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

figure
hist3([data(:,1),data(:,3)],{RANGE,RANGE})
title('Openness vs Extraversion')
xlabel('O')
ylabel('E')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

figure
hist3([data(:,1),data(:,4)],{RANGE,RANGE})
title('Openness vs Agreeableness')
xlabel('O')
ylabel('A')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

figure
hist3([data(:,1),data(:,5)],{RANGE,RANGE})
title('Openness vs Neruoticism')
xlabel('O')
ylabel('N')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

figure
hist3([data(:,2),data(:,3)],{RANGE,RANGE})
title('Conscientiousness vs Extraversion')
xlabel('C')
ylabel('E')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

figure
hist3([data(:,2),data(:,4)],{RANGE,RANGE})
title('Conscientiousness vs Agreeableness')
xlabel('C')
ylabel('A')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

figure
hist3([data(:,2),data(:,5)],{RANGE,RANGE})
title('Conscientiousness vs Neuroticism')
xlabel('C')
ylabel('N')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

figure
hist3([data(:,3),data(:,4)],{RANGE,RANGE})
title('Extraversion vs Agreeableness')
xlabel('E')
ylabel('A')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

figure
hist3([data(:,3),data(:,5)],{RANGE,RANGE})
title('Extraversion vs Neuroticism')
xlabel('E')
ylabel('N')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

figure
hist3([data(:,4),data(:,5)],{RANGE,RANGE})
title('Agreeableness vs Neuroticism')
xlabel('A')
ylabel('N')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
xlabel('A')
ylabel('N')

%%
% Now to actually look at them in 3D space and see if we can't observe any
% clusterings. For example, here are the plots of openness and
% concientiousness vs. extaversion, agreeableness, and neuroticism.

figure
OCE = plot3(data(:,1),data(:,2),data(:,3),'.');

figure
OCA = plot3(data(:,1),data(:,2),data(:,4),'.');

figure
OCN = plot3(data(:,1),data(:,2),data(:,5),'.');

% @@TODO: Figure out how best to visualize density. Change opacity somehow?

%%
% Run k-means on the data. Try various K and examine for how well they
% group, e.g. by plotting groups with colors.

% @@TODO: Plot average distances (and variance of distances?) up to, like,
% k = 50?
clusterList = {};
clusterCent = {};
clusterSum = {};
clusterDist = {};

maxGroups = 50;

for loop = 1:maxGroups
    [clusterList{loop}, clusterCent{loop}, clusterSum{loop}, clusterDist{loop}] = kmeans(data,loop,'start','cluster','replicates',100,'maxiter',10000);
end
%%
% Put the distances between clusters and the centers of the groups they're
% assigned to in a more easily parsed form
clusterSumMeans = [];
for loop = 1:length(clusterSum)
    clusterSumMeans(loop) = mean(clusterSum{loop});
end
clusterSumTotals = clusterSumMeans.*[1:length(clusterSumMeans)];
clusterSumSSE = [];
for loop = 1:75%length(clusterDist)
    clusterSumSSE(loop) = sum(min(clusterDist{loop},[],2).^2);
end

%%
% Find the elbow of the curve by maximizing distance to the linear
% interpolation between first and last points. Taken from
% https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve

%# get coordinates of all the points
nPoints = length(clusterSumSSE);
allCoord = [1:nPoints;clusterSumSSE]';              %'# SO formatting

%# pull out first point
firstPoint = allCoord(1,:);

%# get vector between first and last point - this is the line
lineVec = allCoord(end,:) - firstPoint;

%# normalize the line vector
lineVecN = lineVec / sqrt(sum(lineVec.^2));

%# find the distance from each point to the line:
%# vector between all points and first point
vecFromFirst = bsxfun(@minus, allCoord, firstPoint);

%# To calculate the distance to the line, we split vecFromFirst into two 
%# components, one that is parallel to the line and one that is perpendicular 
%# Then, we take the norm of the part that is perpendicular to the line and 
%# get the distance.
% 
%# We find the vector parallel to the line by projecting vecFromFirst onto 
%# the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
%# We project vecFromFirst by taking the scalar product of the vector with 
%# the unit vector that points in the direction of the line (this gives us 
%# the length of the projection of vecFromFirst onto the line). If we 
%# multiply the scalar product by the unit vector, we have vecFromFirstParallel
scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1), 2);
vecFromFirstParallel = scalarProduct * lineVecN;
vecToLine = vecFromFirst - vecFromFirstParallel;

%# distance to line is the norm of vecToLine
distToLine = sqrt(sum(vecToLine.^2,2));

%# plot the distance to the line
figure('Name','distance from curve to line'), plot(distToLine)

%# now all you need is to find the maximum
[maxDist,idxOfBestPoint] = max(distToLine);

%# plot
figure, plot(clusterSumSSE)
hold on
plot(allCoord(idxOfBestPoint,1), allCoord(idxOfBestPoint,2), 'or')

% This gives results between 5 (for 25 clusters) to 10 (for 75+ clusters).
% Need extra ways to validate.

%%
% @@TODO: perform the same evaluation as above, using the second derivative

%%
% Plot some of the groups and see how they overlap. The axes with highest
% variance are 1,2,5 (openness, conscientiusness, and neuroticism), so
% let's look at those.
close all

for LISTNUM = 3:9
    figure(LISTNUM)
    for loop = 1:10:length(clusterList{LISTNUM})
        switch clusterList{LISTNUM}(loop)
            case 1
                plot3(data(loop,1),data(loop,2),data(loop,5),'.b'), hold on
            case 2
                plot3(data(loop,1),data(loop,2),data(loop,5),'.g'), hold on
            case 3
                plot3(data(loop,1),data(loop,2),data(loop,5),'.r'), hold on
            case 4
                plot3(data(loop,1),data(loop,2),data(loop,5),'.m'), hold on
            case 5
                plot3(data(loop,1),data(loop,2),data(loop,5),'.y'), hold on
            case 6
                plot3(data(loop,1),data(loop,2),data(loop,5),'.k'), hold on
            case 7
                plot3(data(loop,1),data(loop,2),data(loop,5),'.c'), hold on
            case 8
                plot3(data(loop,1),data(loop,2),data(loop,5),'ob'), hold on
            case 9
                plot3(data(loop,1),data(loop,2),data(loop,5),'or'), hold on
            otherwise
                disp('Wrong input list!')
                break
        end
        % drawnow
    end
    title(strcat('Big 5 Clustering:`',strcat(num2str(LISTNUM),' Groups')))
    xlabel('Openness')
    ylabel('Conscientiousness')
    zlabel('Neuroticism')
end

%%
% Evaluating number of clusters with Bayesian Information Criterion
% https://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r
% claims that "BIC = D + ln(n)*m*k" where 
% m = number of columns
% n = number of datapoints
% k = number of centers
% D = total within-cluster sum of squares
for loop = 1:length(clusterDist)
    BIC1(loop) = sum(min(clusterDist{loop},[],2).^2);
    BIC2(loop) = log(length(data))*size(data,2)*loop;
end
figure
plot(BIC1), hold on
plot(BIC2,'k')
plot(BIC1+BIC2,'r'), hold off
% @@TODO implement this correctly, because it appears to be trash currently