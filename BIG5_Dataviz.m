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
%
% One optional modification: subtracting off the mean of each column. (I
% believe, but have not yet checked, that this will make some of the later
% analysis easier, but haven't yet checked.)

% data = data - repmat([mean(data(:,1)),mean(data(:,2)),mean(data(:,3)),mean(data(:,4)),mean(data(:,5))],size(data,1),1);

%%

figure
subplot(2,3,1)
hist(data(:,1),-20:20)
title('Openness')

subplot(2,3,2)
hist(data(:,2),-20:20)
title('Conscientiousness')

subplot(2,3,3)
hist(data(:,3),-20:20)
title('Extraversion')

subplot(2,3,4)
hist(data(:,4),-20:20)
title('Agreeableness')

subplot(2,3,5)
hist(data(:,5),-20:20)
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
% intersection. Instead, construct a 3D histogram:
 
RANGE = [10:50];
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
for loop = 1:length(clusterList)
    clusterSumMeans(loop) = mean(clusterSum{loop});
end
clusterSumTotals = clusterSumMeans.*[1:length(clusterSumMeans)];
clusterSumSSE = clusterSumTotals.^2;

%%
close all
plot(clusterSumMeans,'*-'), hold on
plot(clusterSumTotals,'*-r')
figure
plot(clusterSumSSE,'x-k')

%%
% Visualize the various centroids


%%
% Plot some of the groups and see how they overlap. The axes with highest
% variance are 1,2,5 (openness, conscientiusness, and neuroticism), so
% let's look at those.
close all

for LISTNUM = 3:5
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
    end
    title(strcat('Big 5 Clustering:`',strcat(num2str(LISTNUM),' Groups')))
    xlabel('Openness')
    ylabel('Conscientiousness')
    zlabel('Neuroticism')
end

%%
% Evaluating number of clusters with Bayesian Information Criterion
% https://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r
% claims that "BIC = D + log(n)*m*k" where 
% m = number of columns
% n = number of datapoints
% k = number of centers
% D = total within-cluster sum of squares
for loop = 1:length(clusterDist)
    BIC(loop)= sum(min(clusterDist{loop}).^2) + log10(length(data))*size(data,2)*loop;
end
plot(BIC)

% Currently this is monotone increasing for all clusters 1-90, so either
% there are no real clusters (BIC minimized at k=1) or I've done it wrong.
% Probably the latter.