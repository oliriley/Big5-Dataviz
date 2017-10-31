% Visualizing data from a set of >19k Big Five Inventory responses from an
% anonymized dataset:
% https://openpsychometrics.org/_rawdata/BIG5
close all
clear all

load BIG5_Data.mat
% Each row of the matrix is five numbers representing the sums of a single 
% respondant's answers to the questions of a 50-item Big 5 Inventory. Each 
% of five traits were assessed by 10 questions on a scale from 1-5, so
% values in the matrix will range from 10 (a respondant giving an answer of
% '1' to all questions relevant to a given trait) to 50 (giving an answer
% of '5' to all questions for a trait.
%
% figure
% subplot(2,2,1)
% plot(data(:,1))
% subplot(2,2,2)
% plot(data(:,2))
% subplot(2,2,3)
% plot(data(:,3))
% subplot(2,2,4)
% plot(data(:,4))
% 
% Examining the data shows one error response with 0 across the board.
% There are also a number of responses in which a respondant answered '1'
% or '5' for ALL questions, which strongly suggests they were not taking
% the questionnaire seriously. These have also been removed using:
%
% loop = 1;
% while loop <= length(data(:,1))
%     if (sum(data(loop,:)) <= 50 || (sum(data(loop,:)) >= 250))
%         data = [data(1:loop-1,:) ; data(loop+1:end,:)];
%         loop = loop-1;
%     end
%     loop = loop +1;
% end
% 
% Some anomalous answers remain, but it's good enough for a first pass. 
%
% One optional modification: subtracting off the mean of each column. (I
% believe, but have not yet checked, that this will make some of the later
% analysis easier. Let's find out!)

% data = data - repmat([mean(data(:,1)),mean(data(:,2)),mean(data(:,3)),mean(data(:,4)),mean(data(:,5))],size(data,1),1);

% figure
% subplot(2,3,1)
% hist(data(:,1),-20:20)
% title('Openness')
% 
% subplot(2,3,2)
% hist(data(:,2),-20:20)
% title('Conscientiousness')
% 
% subplot(2,3,3)
% hist(data(:,3),-20:20)
% title('Extraversion')
% 
% subplot(2,3,4)
% hist(data(:,4),-20:20)
% title('Agreeableness')
% 
% subplot(2,3,5)
% hist(data(:,5),-20:20)
% title('Neuroticism')

% Visualizing the histograms of each trait column individually produces
% Gaussian curves, suggesting the traits are normally distributed.
%
% figure
% subplot(3,4,1)
% plot(data(:,1),data(:,2),'.')
% xlabel('O')
% ylabel('C')
% 
% subplot(3,4,2)
% plot(data(:,1),data(:,3),'.')
% xlabel('O')
% ylabel('E')
% 
% subplot(3,4,3)
% plot(data(:,1),data(:,4),'.')
% xlabel('O')
% ylabel('A')
% 
% subplot(3,4,4)
% plot(data(:,1),data(:,5),'.')
% xlabel('O')
% ylabel('N')
% 
% subplot(3,4,5)
% plot(data(:,2),data(:,3),'.')
% xlabel('C')
% ylabel('E')
% 
% subplot(3,4,6)
% plot(data(:,2),data(:,4),'.')
% xlabel('C')
% ylabel('A')
% 
% subplot(3,4,7)
% plot(data(:,2),data(:,5),'.')
% xlabel('C')
% ylabel('N')
% 
% subplot(3,4,8)
% plot(data(:,3),data(:,4),'.')
% xlabel('E')
% ylabel('A')
% 
% subplot(3,4,9)
% plot(data(:,3),data(:,5),'.')
% xlabel('E')
% ylabel('N')
% 
% subplot(3,4,10)
% plot(data(:,4),data(:,5),'.')
% xlabel('A')
% ylabel('N')
%
% Simply plotting each trait against the others is not very illustrative,
% because there is no many indication of how many datapoints lie at each
% intersection. Instead, construct a 3D histogram:
% 
% figure
% hist3([data(:,1),data(:,2)],{-20:20,-20:20})
% title('Openness vs Concientiousness')
% xlabel('O')
% ylabel('C')
% set(gcf,'renderer','opengl')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% 
% figure
% hist3([data(:,1),data(:,3)],{-20:20,-20:20})
% title('Openness vs Extraversion')
% xlabel('O')
% ylabel('E')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% 
% figure
% hist3([data(:,1),data(:,4)],{-20:20,-20:20})
% title('Openness vs Agreeableness')
% xlabel('O')
% ylabel('A')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% 
% figure
% hist3([data(:,1),data(:,5)],{-20:20,-20:20})
% title('Openness vs Neruoticism')
% xlabel('O')
% ylabel('N')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% 
% figure
% hist3([data(:,2),data(:,3)],{-20:20,-20:20})
% title('Conscientiousness vs Extraversion')
% xlabel('C')
% ylabel('E')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% 
% figure
% hist3([data(:,2),data(:,4)],{-20:20,-20:20})
% title('Conscientiousness vs Agreeableness')
% xlabel('C')
% ylabel('A')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% 
% figure
% hist3([data(:,2),data(:,5)],{-20:20,-20:20})
% title('Conscientiousness vs Neuroticism')
% xlabel('C')
% ylabel('N')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% 
% figure
% hist3([data(:,3),data(:,4)],{-20:20,-20:20})
% title('Extraversion vs Agreeableness')
% xlabel('E')
% ylabel('A')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% 
% figure
% hist3([data(:,3),data(:,5)],{-20:20,-20:20})
% title('Extraversion vs Neuroticism')
% xlabel('E')
% ylabel('N')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% 
% figure
% hist3([data(:,4),data(:,5)],{-20:20,-20:20})
% title('Agreeableness vs Neuroticism')
% xlabel('A')
% ylabel('N')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% xlabel('A')
% ylabel('N')

% Now to actually look at them in 3D space and see if we can't observe any
% clusterings. For example, here are the plots of openness and
% concientiousness vs. extaversion, agreeableness, and neuroticism.

% figure
% OCE = plot3(data(:,1),data(:,2),data(:,3),'.');
% 
% figure
% OCA = plot3(data(:,1),data(:,2),data(:,4),'.');
% 
% figure
% OCN = plot3(data(:,1),data(:,2),data(:,5),'.');

% @@TODO: Figure out how to visualize density

% Run k-means on the data. Try various K and examine for how well they
% group, e.g. by plotting groups with colors.

% [clusterList9, clusterCent9] = kmeans(data,9,'start','cluster','replicates',10,'maxiter',1000);
% figure
% for loop = 1:length(clusterList9)
%     if clusterList9(loop) == 1
%         plot3(data(loop,1),data(loop,2),data(loop,3),'.b'), hold on
%     end
%     if clusterList9(loop) == 3
%         plot3(data(loop,1),data(loop,2),data(loop,3),'.g'), hold on
%     end
%     if clusterList9(loop) == 9
%         plot3(data(loop,1),data(loop,2),data(loop,3),'.r'), hold on
%     end
% end
% plot3(clusterCent9(:,1),clusterCent9(:,2),clusterCent9(:,3),'.')

[clusterList16, clusterCent16] = kmeans(data,16,'start','cluster','replicates',10,'maxiter',1000);
figure
for loop = 1:length(clusterList16)
    if clusterList16(loop) == 1
        plot3(data(loop,1),data(loop,2),data(loop,3),'.b'), hold on
    end
    if clusterList16(loop) == 3
        plot3(data(loop,1),data(loop,2),data(loop,3),'.g'), hold on
    end
    if clusterList16(loop) == 9
        plot3(data(loop,1),data(loop,2),data(loop,3),'.r'), hold on
    end
end
plot3(clusterCent16(:,1),clusterCent16(:,2),clusterCent16(:,3),'.')

% @@TODO: Plot average distances (and variance of distances?) across, like,
% k = 2:100? 100 is suggested by rule of thumb k ~ sqrt(n/2)