function ScoreMatrix = DCF(split, k, featureRank, networkRank, lambda, alpha)
    %%%% INPUT PARAMS:
    %%% split       : One of the k-fold splits to be hidden.
    %%% k           : Rank of the parameter matrix Z.
    %%% featureRank : Number of compressed dimensions (microarray gene expression data).
    %%% networkRank : Number of latent features from similarity networks.
    %%% lambda      : Regularization parameter l in the objective.
    %%% alpha       : the parameter r for PU formulation

    %%%% OUTPUT:
    %%% ScoreMatrix : Matrix of size # genes x # diseases.
    
	%%% Data files are provided i n the same directory.
    load ('genes_phenes.mat');
    load ('GeneFeatures.mat');
    load ('clinicalfeatures_tfidf.mat');
    
	GenePheneTraining = GenePhene{1};%% take one of GenePhene as train set
	if (split > 0) %% else use all data to make predictions.
	    load ('splits_uniform.mat');
	    numPhenes = size(splits{split},2); 
	    GenePheneTraining = GenePheneTraining(:,1:numPhenes) - splits{split};%% 
	end
	Features = [];
	ColFeatures = [];
        
	if (featureRank == Inf)  
        Features = GeneFeatures; 
        ColFeatures = F;
	elseif (featureRank > 0) 
        display ('Reducing feature dimensionality..');
        %% Reducing dimensionality of microarray
        GeneFeatures = normalization(GeneFeatures);
		mysdae(GeneFeatures, featureRank, 0); 
         load('sdae999.mat');
         Features = H;
         clear H;
       %% Reducing dimensionality of OMIM pages
        [U,S] = svds(sparse(F(1:numPhenes,:)), featureRank);  		
        ColFeatures = U;
  	end 

 	if (networkRank > 0) 
 		display ('Reducing gene network dimensionality..');
        %% Reducing dimensionality of gene-gene interaction network
 		mysdae(full(GeneGene_Hs), networkRank, 1); 
 		 load('sdae999.mat');   
         Features = [Features H(1:numGenes,:)];
         clear H
 	    %% Reducing dimensionality of orthologous phenotypes
  		GP_sp = []; 
      		for sp=2:numel(GenePhene) 
          		GP = GenePhene{sp};  
          		for i=1:size(GP,2)   
              		if(sum(GP(:,i)) > 0)               		
                        GP(:,i) = GP(:,i)/norm(GP(:,i)); 
                    end
                end
                GP_sp = [GP_sp GP];
            end ;
 		mysdae(full(GP_sp), networkRank, 1);
         load('sdae999.mat');
         Features = [Features H]; 
         clear H
        %% Reducing dimensionality of  disease similarities network
        [U,S] = svds(PhenotypeSimilaritiesLog(1:numPhenes,1:numPhenes), networkRank);
  		ColFeatures = [ColFeatures U];       
 		clear U S
        end
 
	for i=1:size(Features,2)
		if (norm(Features(:,i)) > 0)
			Features(:,i) = Features(:,i)/norm(Features(:,i));
		end
	end
	for i=1:size(ColFeatures,2)
		if (norm(ColFeatures(:,i))>0)
			ColFeatures(:,i) = ColFeatures(:,i)/norm(ColFeatures(:,i));
		end
	end
        %% If no features are supplied (i.e. networkRank & featureRank are 0), use default matrix completion.
	if (isempty(Features))
		Features = eye(numGenes); 
	end
	if (isempty(ColFeatures))
		ColFeatures = eye(numPhenes);
	end
	fprintf('Using %d row, %d col features.\n', size(Features,2), size(ColFeatures,2));
  
    [W H rmse walltime] = imf_train(sparse(double(GenePheneTraining)), sparse(double(Features)), sparse(double(ColFeatures)),  ['-n 16 -t 15 -T 5 -g 20 -p 3 -a 0 -s 11 -k ', num2str(k), ' -l ', num2str(lambda),  ' -r ', num2str(alpha)]);
    ScoreMatrix = Features * W *H' * ColFeatures';
 end
