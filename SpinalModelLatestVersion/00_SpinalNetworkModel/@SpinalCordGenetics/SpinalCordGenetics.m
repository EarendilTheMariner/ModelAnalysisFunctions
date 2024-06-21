classdef SpinalCordGenetics 
    methods (Static) 
        function Genetics = GetMouseGenetics() 
        % if(~exist('ConProbCentral'))
        % %% Directional Probability
        %     [ConProbCentral,ModesC] = ComputeModeSignificance(Guidance(GCS&GCT,:),categorical([Central]),categorical([Central]),'Perm',1,'P-Val',0.05,'Mode','Pair'); %Pair
        %     [ConProbSensory,ModesS] = ComputeModeSignificance(Guidance(GCS&GST,:),categorical([Central]),categorical([Sensory]),'Perm',1,'P-Val',0.05,'Mode','Pair'); %Pair
        % end
        % %%
        % if(~exist("GeneConProb"))
        %     [GeneConProb,CP] = ComputeSynapseProba(ConProbCentral,Res3D,'Discrim','Class');
        %     [GeneConProbS,CPS] = ComputeSynapseProba(ConProbSensory,Res3D,'Discrim','Class');
        % end
            Genetics = nan;
        end
    end
end
