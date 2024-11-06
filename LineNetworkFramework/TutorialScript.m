%% Quick tutorial for building 1D networks with 5 subpopulations

% Instantiate a network object

MyNetwork = LineNetwork(Length,PopParams);

% Length is the spatial length of the network. Maybe stick to 1000 for now.
% I also specifies the number of neurons, keeping the spatial density of
% neurons fixed.

% PopParams is a 1x5 vector, where each element is between 0 and 1. It
% specifies the how far along the normalized pair-wise distance domain each
% subpopulation projects on average. It determines the spatial projection
% probability distributions for each subpopulation. For example, if
% PopParams = [0.1, 0.3, 0.5, 0.7, 0.9], it means that subpopulation (sp) 1
% projects very caudally, sp2 projects caudally at more intermediate
% distances, sp3 projects exactly locally, sp4 projects rostrally
% intermediate and sp5 projects far rostrally. sp1, sp3 and sp5 are
% excitatory, sp2 and sp4 are inhibitory. 

% You can get an overview of the network summary statistics by plotting:

figure;
SpatialCoupling(MyNetwork.ConnMat,MyNetwork.E,MyNetwork.I,MyNetwork.Coordinates,true);

% This is the "projectome" of the instantiated network.

% Once you've built the network, you can use two methods:

MyNetwork.SimulateLine(t_steps);

% This will simulate the network and store the rates in the MyNetwork.Rates property,
% as a Neuron x Time matrix. t_steps is the simulation duration.  

MyNetwork.GainMod(Pop,Gain);

% This will scale the synapses of the population of interest. If you want
% to double the strength of synapses belonging to subpopulation 3, you
% would put Pop = 3 and Gain = 2. Simple. 

% Start by reading the LineNetwork object class, understand the Properties
% and Methods. I've also added a folder of UtilityFunctions, which might be
% helpful for computing different stats and metrics on your network. It
% also includes a function for building random networks, which you were
% introduced to in my thesis. 
% 
% Feel free to always ask me questions, there are many ways
% this code can be optimized lol, so it will for sure break down if you
% test the argument ranges
