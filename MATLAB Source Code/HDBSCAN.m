
classdef HDBSCAN < handle
    % clusterer = HDBSCAN( X )
    %
    % Calling HDBSCAN creates an instance of the HDBSCAN cluster object.
    %
    % HDBSCAN stands for: Hierarchical Density-Based Spatial Clustering, 
    % with Application with Noise. It is extensively described in:
    %   Campello et al. 2013 and Campello et al. 2015
    %
    % The HDBSCAN cluster object contains methods for training a hierarchical
    % clustering model based on the input data matrix X. Clustering is
    % performed by iteratively removing links between a
    % graph-representation of the original data (based on pairwise
    % distances), and searching for the resulting connected
    % components at each iteration. 
    %
    % This differs from other hierarchical clustering methods,
    % such as single linkage clustering, as clusters with less
    % than a minimum # of points are deemed noise. Additionally, following
    % the creation of the cluster hierarchy, an optimal, flat clustering is
    % performed based on the stability of each cluster. This gives a final 
    % clustering that can be performed at varying heights for each branch
    % of the cluster tree, which differs than DBSCAN which produces only a
    % single, horizontal cut through the tree.
    %
    % After training a model, one can also predict cluster membership for new
    % points not originally used in the model creation. Note that this gives an 
    % "approximate membership", as new points may have changed the model
    % hierarchy if used in the training procedure.
    %
    % Properties:
    % ----------
    %   data            -   the raw data used for model creation
    %
    %   nPoints         -   the number of rows in the matrix "data"
    %
    %   nDims           -   the number of columns in the matrix "data"
    %   
    %   minpts          -   the nearest 'minpts' neighbor used for core distance
    %                       calculation for each point in X. Default = 5
    %   
    %   minclustsize    -   the minimum # of points necessary for a cluster
    %                       to be deemed valid. Default = 5
    %
    %   minClustNum     -   the minimum # of clusters to be realized. Default = 1
    %
    %   outlierThresh   -   a cutoff value between [0,1], where any X(i) with an outlier score 
    %                       (see below) greather than 'outlierThresh' is assigned
    %                       as an outlier (ID = 0). Default = 0.9
    %
    %   kdtree          -   a KD tree based on the data in X if the
    %                       dimensionality of X is <= 10
    %
    %   model           -   a trained hierarchical model. For details, see
    %                       'hdbscan_fit.m'.
    %
    %   labels          -   vector of cluster membership for each point i
    %
    %   bestClusters    -   the optimal clusters discovered from the clusterTree
    %
    %   clusterMap      -   maps each unique ID in labels to the best
    %                       cluster it is associated with
    %
    %   corePoints      -   the most "representative" points of the final
    %                       optimal clusters. These are the densest points
    %                       of any of the clusters, and can be used for
    %                       predicting new data cluster membership
    %
    %   coreLambda      -   the lambda value associated with the core
    %                       points for each best cluster
    %
    %   score           -   the outlier score for each point i
    %
    %   dCore           -   the core distances of the points in X, given
    %                       the specified 'minpts'
    %
    %   P               -   probability of point i belonging to the cluster
    %                       indicated by labels(i)
    %
    % Methods:
    % -------
    %   fit_model       -   fits a hierarchical model to the data in X
    %
    %   predict         -   predicts new data based on the trained model
    %
    %   get_best_clusters - finds the optimal flat clustering from the
    %                       full cluster_tree hierarchy
    %
    %   get_membership  -   assigns a label and probability to each point in X 
    %                       based on the optimal flat clustering.
    %
    %   plot_tree       -   plots the cluster hierarchy, and indicates
    %                       which clusters were kept in the final clustering
    %
    %   plot_clusters   -   plots the first 3 (or 2, if 2D) columns of self.data,
    %                       color coded by the cluster labels of the data points
    %
    %   run_hdbscan     -   convenience function that fits a full
    %                       hierarchical model, finds optimal clusters, and
    %                       assigns labels to data points
    %
    %       * see individual methods for more details on inputs/outputs
    %
    % Written by Jordan Sorokin
    % 10/15/2017
    
    properties
        nPoints
        nDims
        model
        kdtree
        data
        minpts = 5;
        minclustsize = 5;
        minClustNum = 1;
        outlierThresh = 0.9;
        bestClusters
        clusterMap
        corePoints
        coreLambda
        labels
        score
        P
    end
    
    methods
        
        function self = HDBSCAN( X )
            % creates an instance of the HDBSCAN object
            
            self.data = X;
            self.nPoints = size( X,1 );
            self.nDims = size( X,2 );
        end
            
        
        function fit_model( self,varargin )
            % fit_model( self,(dEps,verbose) )
            %
            % fits a full hierarchical cluster model to the data stored in 
            % self.data. Uses "self.minpts" and "self.minclustsize" for
            % training the model. 
            %
            % Inputs:
            %   self - an instance of the HDSBCAN object
            %
            %   dEps - a scalar that specifies the number of iterations to
            %          do, as:
            %           
            %               nIters = iterations(1:dEps:end)
            %
            %           Larger values of 'dEps' results in faster model
            %           training, at the risk of more approximate cluster
            %           identification (default = 1)
            %   
            %   verbose - logical. prints clustering information if true
            %
            % Outputs:
            %   self.model
            
            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                dEps = round( varargin{1} ) * sign( varargin{1} ); % ensures it's a positive integer
            else
                dEps = 1;
            end
            
            if nargin > 2 && ~isempty( varargin{2} )
                verbose = varargin{2};
            else
                verbose = true;
            end
            
            % remove previous cluster-based post processing
            self.bestClusters = [];
            self.corePoints = [];
            self.coreLambda = [];
            self.P = [];
            self.score = [];
            self.labels = [];
            
            % report cluster params if verbose = true
            if verbose
                fprintf( 'Training cluster hierarchy...\n' );
                fprintf( '\tData matrix size:\n' );
                fprintf( '\t\t%i points x %i dimensions\n\n',self.nPoints,self.nDims );
                fprintf( '\tMin # neighbors: %i\n',self.minpts );
                fprintf( '\tMin cluster size: %i\n',self.minclustsize );
                fprintf( '\tMin # of clusters: %i\n',self.minClustNum );
                fprintf( '\tSkipping every %i iteration\n\n',dEps-1 );
                start = clock;
            end
                
            % fit the hierarchical cluster tree
            self.model = hdbscan_fit( self.data,...
                                'minpts',self.minpts,...
                                'minclustsize',self.minclustsize,...
                                'minClustNum',self.minClustNum,...
                                'dEps',dEps );
            
            % report time to fit the model               
            if verbose
                stop = clock;
                fprintf( 'Training took %0.3f seconds\n',(stop(end-1)*60+stop(end)) - (start(end-1)*60+start(end)) );
            end
        end
        
        
        function get_best_clusters( self )
            % get_best_clusters( self )
            %
            % produces the optimal flat clustering from the hierarchical
            % cluster scheme in self.model by finding the most stable 
            % clusters in a recusive way
            %
            % Outputs:
            %   self.bestClusters
            %   self.corePoints
            %   self.coreLambda
            
            % check if model has been trained
            self.trained_check();
            tree = self.model.clusterTree;
            
            % get the optimal flat clustering
            self.bestClusters = OFCH( tree.stability,tree.parents );
            
            % find maximum lambda and core points for the best clusters
            [self.corePoints,self.coreLambda] = get_core_points( tree.parents,self.bestClusters,full( self.model.lambdaMax ) );
        end  
        
        
        function get_membership( self )
            % get_membership( self )
            %
            % finds the cluster membership of the points in self.data based
            % on the best clusters found in the hierarchy
            %
            % Outputs:
            %   self.labels
            %   self.score
            %   self.P
            
            % check if model has been trained           
            self.trained_check();
            
            % check if we've performed optimal flat clustering
            self.best_cluster_check();

            % compute the outlier scores
            tree = self.model.clusterTree;
            self.score = GLOSH( self.bestClusters,tree.parents,self.model.lastClust,self.model.lambdaNoise );
            
            % compute labels and probability of cluster membership
            [self.labels,self.P] = get_cluster_probability( self.bestClusters,full( self.model.lambdaMax ),self.coreLambda );
            self.clusterMap = unique( self.labels(self.labels>0) )';

            % set labels with outlier scores > outlierThresh = 0
            self.labels( self.score > self.outlierThresh ) = 0;    
            
            % update if any clusters are now all zero
            badclusts = ~ismember( self.clusterMap,unique( self.labels(self.labels>0) ) );
            self.clusterMap( badclusts ) = [];
            self.corePoints = self.corePoints( self.clusterMap );
            self.coreLambda = self.coreLambda( self.clusterMap );
            self.bestClusters = self.bestClusters( self.clusterMap );
        end
        
        
        function run_hdbscan( self,varargin )
            % run_hdbscan( self,(minpts,minclustsize,minClustNum,dEps,outlierThresh,plotResults) )
            %
            % fits a hierarchical model to self.data and finds the best
            % flat clustering scheme. Then assigns labels to each data
            % point in self.data based on the final clusters.
            %
            %   Note: this is just a convenience function to avoid manually
            %         typing the commands to perform these operations
            %
            % Inputs:
            %   minpts - minimum # neighbors for core distances
            %
            %   minclustsize - minimum # points in a cluster to keep the
            %                  cluster
            %
            %   minClustNum - the minimum # of clusters to be discovered. Default = 1
            %
            %   dEps - iterator (skips edge weight iteratios every "dEps"
            %          times)
            %
            %   outlierThresh - threshold between [0,1] for outlier scores
            %
            %   plotResults - logical to plot the cluster results or not
            %
            % Outputs:
            %   self.model
            %   self.corePoints
            %   self.coreLabels
            %   self.bestClusters
            %   self.labels
            %   self.P
            %   self.score
            
            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                self.minpts = varargin{1};
            end
            if nargin > 2 && ~isempty( varargin{2} )
                self.minclustsize = varargin{2};
            end
            if nargin > 3 && ~isempty( varargin{3} )
                self.minClustNum = varargin{3};
            else
                self.minClustNum = 1;
            end
            if nargin > 4 && ~isempty( varargin{4} )
                dEps = varargin{4};
            else
                dEps = 1;
            end
            if nargin > 5 && ~isempty( varargin{5} )
                self.outlierThresh = varargin{5};
            end
            if nargin > 6 && ~isempty( varargin{6} )
                plotResults = varargin{6};
            else
                plotResults = false;
            end
            
            % fit the hierarchical model
            self.fit_model( dEps );
            
            % extract best clusters
            self.get_best_clusters();
            
            % assign labels
            self.get_membership();
            
            % visualize the results
            if plotResults
                figure;
                self.plot_tree();
                
                figure;
                self.plot_clusters();
                set( gcf,'color','k' )
            end
        end
        
        
        function update_hierarchy( self,newLabels )
            % update_hierarchy( self,newLabels )
            %
            % updates the cluster hierarchy and the lambda values associated
            % with the clusters based on any new label vector. This allows
            % one to manually alter the clusters while maintaining a
            % probabilistic model that can be used to predict new points
            %
            % Inputs:
            %   newLabels - self.nPoints x 1 vector of new labels
            %
            % Outputs:
            %   updates all properties pertaining to the cluster hierarchy
            %   in "self"
            
            % check for model / best clusters
            self.trained_check();
            self.best_cluster_check();
            if isrow( newLabels )
                newLabels = newLabels';
            end
            
            % get the necessary variables that will be updated
            lambdaMax = full( self.model.lambdaMax );
            lambdaNoise = self.model.lambdaNoise;
            bestClusts = self.bestClusters;
            map = self.clusterMap;
            parents = self.model.clusterTree.parents;
            clusters = self.model.clusterTree.clusters;
            minLambda = self.model.clusterTree.lambdaMin;
            stability = self.model.clusterTree.stability;
            nClusts = clusters(end);
            newClusters = [];
            
            % find changed labels
            oldLabels = self.labels;
            changedPts = (oldLabels ~= newLabels);
            changedLabels = unique( newLabels(changedPts) )';
            changedLabels(changedLabels == 0) = [];

            % loop over the changed clusters, and update the model 
            % depending on whether the new cluster is a result 
            % of a split or merge
            for k = changedLabels
                pts = changedPts & (newLabels == k); % intersection{ C_k, C_i }
                prevID = unique( oldLabels(pts) );
                prevClust = bestClusts( ismembc( map(~ismembc(bestClusts,newClusters)),prevID ) );
                thisClust = bestClusts( map == k );
                
                switch numel( prevClust )
                    case 0 % previous points were just noise
                        
                        % update lambdas of points by using the lambdas at
                        % which they become noise
                        nClusts = nClusts + 1;
                        newLambda = lambdaNoise(pts);
                        lambdaMax(:,nClusts) = 0;
                        lambdaMax(pts,nClusts) = newLambda;
                        
                        % update cluster stability/lambda by taking the
                        % mean of the parent clusters
                        clusters(nClusts) = nClusts;
                        parents(nClusts) = max( round( mean( self.model.lastClust(pts) ) ),1 ); % this is a hack
                        minLambda(nClusts) = minLambda(parents(nClusts));
                        stability(nClusts) = stability(parents(nClusts));
                        bestClusts(end+1) = nClusts;
                        map(end+1) = k;
                        newClusters(end+1) = nClusts;
                                  
                    case 1 && ~any( thisClust ) % split / manual new cluster
                        
                        % update the lambdas of the new clusters by simply
                        % moving lambdas associated with appropraite points
                        nClusts = nClusts + 1;
                        newLambda = lambdaMax(pts,prevClust);
                        lambdaMax(pts,prevClust) = 0;
                        lambdaMax(:,nClusts) = 0;
                        lambdaMax(pts,nClusts) = newLambda;
                        
                        % update cluster minimum lambda by just copying
                        % from previous one
                        minLambda(nClusts) = minLambda(prevClust); % just take a copy
                        stability(nClusts) = stability(prevClust); % ditto
                        
                        % update the clusters and parents vectors
                        clusters(nClusts) = nClusts;
                        parents(nClusts) = prevClust;
                        bestClusts(end+1) = nClusts;
                        map(end+1) = k;
                        newClusters(end+1) = nClusts;

                    otherwise % merge
                        
                        % check if the current cluster is a new cluster,
                        % resulting from a merge that produced an entirely
                        % new ID, rather than merging with a previous ID
                        if isempty( thisClust )
                            nClusts = nClusts + 1;
                            bestClusts(end+1) = nClusts;
                            map(end+1) = k;
                            newClusters(end+1) = nClusts;
                            thisClust = nClusts;
                            lambdaMax(:,thisClust) = 0;
                            minLambda(thisClust) = 0;
                            stability(thisClust) = 0;
                        end
                        
                        for j = 1:numel( prevID )
                            
                            % update the lambdas
                            oldpts = (oldLabels == prevID(j));
                            ptFrac = nnz( oldpts ) / nnz( newLabels==k );
                            switch prevID(j) 
                                case 0 % just noise
                                    %oldCluster = mode( self.model.lastClust(oldpts) );
                                    oldLambda = lambdaNoise(oldpts);
                                    
                                otherwise % prev cluster used
                                    oldCluster = bestClusts(map == prevID(j) & ~ismember( bestClusts,newClusters ));
                                    oldLambda = lambdaMax(oldpts,oldCluster);
                                    lambdaMax(oldpts,oldCluster) = 0;
                                    
                                    % update the minLambda and stability vectors
                                    % by taking a weighted average, determined by
                                    % the fraction of points merged from cluster j
                                    minLambda(thisClust) = (minLambda(thisClust) + ptFrac*minLambda(oldCluster)) / 2;
                                    stability(thisClust) = (stability(thisClust) + ptFrac*stability(oldCluster)) / 2;
                            end
                            
                            lambdaMax(oldpts,thisClust) = oldLambda;
                        end 
                end
            end
                        
            % eliminate old clusters that are now just noise or eliminated
            newClusts = unique( newLabels(newLabels > 0) );
            badClusts = ~ismembc( map,newClusts );
            bestClusts(badClusts) = [];
            map(badClusts) = [];
            
            % now find the core points and core lambda for the clusters
            [self.corePoints,self.coreLambda] = get_core_points( parents,bestClusts,lambdaMax );

            % store the updated parameters
            for i = 1:numel( map )
                newLabels(newLabels==map(i)) = i;
                map(i) = i;
            end
            self.bestClusters = bestClusts;
            self.clusterMap = map;
            self.model.lambdaMax = sparse( lambdaMax );
            self.model.clusterTree.clusters = clusters;
            self.model.clusterTree.parents = parents;
            self.model.clusterTree.stability = stability;
            self.model.clusterTree.lambdaMin = minLambda;
            
            % update the labels
            self.labels = newLabels;
        end
        
        
        function [newLabels,newProb,outliers] = predict( self,newPoints )
            % [newLabels,newProb,outliers] = predict( self,newPoints,alpha )
            %
            % predicts cluster membership to new points given the trained
            % hierarchical cluster model.
            %
            % For each point i, prediction is performed as follows:
            %   (a) set D(i,j) = euclidean distance of the jth nearest
            %       neighbor with label > 0, for j in [1, self.minpts*2]
            %
            %   (b) set R(i) = nearest self.minpts mutual-reachability 
            %       distance among the self.minpts*2 nearest neighbors.
            %       
            %   (c) assign label(i) = label of the nearest mutual-reachable
            %       neighbor of point i
            %   
            %   (d) set L(i) = lambda value for point i as 1 / R(i)
            %
            %   (e) P(i) = L(i) / L_max; L_max = maximum lambda of the
            %       cluster assigned to point i
            %
            %   (f) flag outlier(i) IF 1-P(i) > self.outlierThresh
            
            % check trained and best cluster assignment
            self.trained_check()
            self.best_cluster_check();
            
            % check sizes of data matrices
            [n,m] = size( newPoints );
            assert( m == self.nDims,'new data must have same # of columns as trainng data' );
            
            % create a kdtree object for finding nearest neighbors
            if isempty( self.kdtree )
                self.create_kdtree();
            end

            % for each point, find the nearest core points 
            uID = unique( self.labels(self.labels>0) );
            nID = numel( uID );
            D = zeros( n,nID ); 
            for i = 1:nID
                points = self.data( self.corePoints{i},: );
                d = compute_pairwise_dist( newPoints,points );
                D(:,i) = min( d,[],2 );
            end
            
            % convert the mutual reaches to lambda values
            [newLambda,newLabels] = min( D,[],2 );
            newLambda = 1./newLambda;
            
            % now that we have the lambda values, we can check if any of
            % the new points are outliers, by comparing their lambda values
            % with the minimum lambda value of the clusters they are
            % assigned to. This relates to the largest "weight" in the
            % original hierarchical tree that a point can have 
            % while still being associated with its particular cluster
            uniqueLabels = unique( newLabels(newLabels>0) )';
            newProb = zeros( size( newLabels ) );
            lambdaCore = self.coreLambda;
            map = self.clusterMap;

            % compare the lambda values to the max lambda of this
            % cluster (the core points) to get the probability of
            % belonging to this cluster
            for k = uniqueLabels
                thesePts = (newLabels == k);
                newProb(thesePts) = newLambda(thesePts) ./ max( lambdaCore(map == k),newLambda(thesePts) );
            end
            
            % outlier if 1 - probability is > outlier threshold
            outliers = find( (1-newProb) > self.outlierThresh );
        end

        
        function plot_tree( self )
            % plot_tree( self )
            % 
            % plots the cluster hierarchy tree stored in self.model
            
            % check if trained
            trained_check( self )
            
            % create the plot and change plot style
            [~,h] = plot_cluster_tree( self.model.clusterTree );
            nclusts = length( self.model.clusterTree.clusters );
            h.MarkerSize = 4;
            h.NodeColor = repmat( [0 0 0],nclusts,1 );
            h.LineStyle = '--';
            h.EdgeColor = 'k';
            h.NodeLabel = repmat( {''},1,nclusts );
            set( gca,'tickdir','out','box','off','XTick',[],'XTickLabel',[] );
            title( 'Condensed cluster tree' );

            % highlight kept clusters
            if ~isempty( self.bestClusters )
                h.NodeColor(self.bestClusters,:) = repmat( [1 0 0],length( self.bestClusters ),1 );
                h.NodeLabel(self.bestClusters) = strsplit( num2str( self.bestClusters ),' ' );
            end
        end
        
        
        function h = plot_clusters( self,varargin )
            % h = plot_clusters( self,(dims) )
            %
            % plots the clusters, color-coded by the labels,
            % defaulting to the first 3 columns of self.data
            %
            % Inputs:
            %   (dims) - up to 3 dimensions (columns) of self.data to plot.
            %            Must specify self.nDims different dims to plot
            %
            % Outputs:
            %   h - handle to scatter plot
            
            if nargin > 1 && ~isempty( varargin{1} )
                dims = varargin{1};
                dims = dims(1:min( self.nDims,3 ));
            else
                dims = 1:self.nDims;
            end
            
            % scatter plots
            if self.nDims >= 3
                h = scatter3( self.data(:,dims(1)),self.data(:,dims(2)),self.data(:,dims(3)),'.' );
            else
                h = scatter( self.data(:,dims(1)),self.data(:,dims(2)),'.' );
            end
            
            % change colors according to self.labels
            if ~isempty( self.labels )
                h.CData = self.labels;
                colormap( self.cluster_colors );
            end
            
            % change appearance
            title( 'Clustered data','color','w' );
            set( h.Parent,'tickdir','out','box','off','color','k','xcolor','w','ycolor','w' );
        end

    end % public methods
    
    
    %% private methods
    methods(Access=private)
        
        function trained_check( self )
            % returns an error if self.trained is false
            
            assert( ~isempty( self.model ),'Must train hierarchical model first!' );
        end
        
        function best_cluster_check( self )
            % returns an error if self.bestClusters is empty
            
            assert( ~isempty( self.bestClusters ),'No optimal flat clusters found.' );
        end
        
        function colors = cluster_colors( self )

             plotColor = [ 
                 [.65, .65, .65];...   % light gray         (0)
                 [0.1, 0.74, 0.95];...  % deep sky-blue     (1)
                 [0.95, 0.88, 0.05];... % gold/yellow       (2)
                 [0.80, 0.05, 0.78];... % magenta           (3)
                 [0.3, 0.8, 0.20];...   % lime green        (4)
                 [0.95, 0.1, 0.1];...   % crimson red       (5)   
                 [0.64, 0.18, 0.93];... % blue-violet       (6)
                 [0.88, 0.56, 0];...    % orange            (7)
                 [0.4, 1.0, 0.7];...    % aquamarine        (8)
                 [0.95, 0.88, 0.7];...  % salmon-yellow     (9)
                 [0, 0.2, 1];...        % blue              (10)
                 [1, 0.41, 0.7];...     % hot pink          (11)
                 [0.5, 1, 0];...        % chartreuse        (12)
                 [0.6, 0.39, 0.8];...   % amtheyist         (13)
                 [0.82, 0.36, 0.36,];...% indian red        (14)
                 [0.53, 0.8, 0.98];...  % light sky blue    (15)
                 [0, 0.6, 0.1];...      % forest green      (16)
                 [0.65, 0.95, 0.5];...  % light green       (17)
                 [0.85, 0.6, 0.88];...  % light purple      (18)
                 [0.90, 0.7, 0.7];...   % light red         (19)
                 [0.2, 0.2, 0.6];...    % dark blue         (20)
                ];

            repeats = max( 1,ceil( max( self.labels )/19 ) );
            colors = [plotColor(1,:);repmat( plotColor(2:end,:),repeats,1 )];
        end
        
        function create_kdtree( self )
            % creates a kdtree object based on self.data. Used for fast
            % nearest-neighbor queries for new points
            
            self.kdtree = createns( self.data(self.labels>0,:),'nsmethod','kdtree' );
        end

    end

end 


function score = GLOSH( clusters,parent,lastClust,lambdaNoise )
    % Global-Local Outlier Score from Hierarchy 
    %
    % score = GLOSH( delta,parent,lastClust,lambdaNoise )
    %
    % Computes the outlier score for each point 1:n given a 
    % cluster indicator vector "delta" and minimum edge weight 
    % at which each point i existed in the hierarchy
    %
    % Inputs
    %   clusters - 1 x K vector of cluster IDs
    %
    %   parent - 1 x K vector of parent clust for each cluster i
    %
    %   lastClust - 1 x n vector specifying the last cluster i that each
    %               point belonged to before becoming noise
    %
    %   lambdaNoise - 1 x n vector specifying the lambda value beyond which
    %                 point i becomes noise
    %
    % Output
    %   score - 1 x n vector of outlier scores, computed as:
    %   
    %               score(i) = 1 - refWeight / lastWeight(i)
    %
    %           where refWeight = min( lastWeight )
    % 
    % Written by Jordan Sorokin, 10/5/17
    
    score = zeros( size( lastClust ) );
    hasProcessed = false( 1,max( clusters ) );
    
    for clust = clusters
        
        % check if we've processed this cluster
        if hasProcessed(clust)
            continue
        end
        
        % perform depth-first search from this node to get children
        children = get_all_subnodes( clust,parent );

        % get epsilon when this cluster or any subclusters completely
        % disappear (when all X(i) == 0 for this branch)
        pts = ismember( lastClust,[clust,children] );
        refLambda = max( lambdaNoise(pts) );

        % compute the Global-Local Outlier Scores
        score(pts) = 1 - (lambdaNoise(pts) ./ refLambda); % [EQ 8]
        
        hasProcessed([clust,children]) = true;
    end
end

function [delta,S_hat] = OFCH( S,parent )
    % Optimal Flat Clustering from Hierarchy
    %
    % [delta,S_hat] = OFCH( S,parentClust )
    %
    % Computes the optimal cluster scheme from a cluster hierarchy
    % cluster stability vector. 
    %
    % Inputs
    %   S - 1 x K vector of cluster stabilities, where K = # clusters
    %   
    %   parent - 1 x K vector of parent clusters of current cluster.
    %            i.e. parent(i) = parent of ith cluster
    %
    % Outputs
    %   delta - 1 x m < K vector of optimal clusters
    %
    %   S_hat - 1 x m vector of modified cluster stabilities
    %           
    %               S_hat(i) = max{ S(i), sum( S_hat(c) ) }
    %           
    %           for c in the set of all children of cluster i
    % 
    % Written by Jordan Sorokin, 10/5/2017

    % get the parent-child node pair and preallocate indicator vectors
    uniqueParent = fliplr( unique( parent(parent>0) ) );
    S_hat = S;
    delta = true( 1,numel( S ) );

    for clust = uniqueParent

        % sum stability of children of current parent
        children = (parent == clust);
        S_child = sum( S_hat(children) );

        % compare with stability of parent and store the larger
        [S_hat(clust),ind] = max( [S(clust),S_child] ); % [EQ 5] 

        % update our indicator vector
        if ind == 1
            % this parent node is stable, so set children to 0
            allChildren = get_all_subnodes( clust,parent );
            delta(allChildren) = false;
        else
            % children more stable, so set current clust to 0
            delta(clust) = false;
        end
    end
    
    delta = find( delta );
end

function [d dt pred] = bfs(A,u,target)
% BFS Compute breadth first search distances, times, and tree for a graph
%
% [d dt pred] = bfs(A,u) returns the distance (d) and the discover time
% (dt) for each vertex in the graph in a breadth first search 
% starting from vertex u.
%   d = dt(i) = -1 if vertex i is not reachable from u
% pred is the predecessor array.  pred(i) = 0 if vertex (i)  
% is in a component not reachable from u and i != u.
%
% [...] = bfs(A,u,v) stops the bfs when it hits the vertex v
%
% Example:
%   load_gaimc_graph('bfs_example.mat') % use the dfs example from Boost
%   d = bfs(A,1)
%
% See also DFS

% David F. Gleich
% Copyright, Stanford University, 2008-20098

% History
% 2008-04-13: Initial coding

if ~exist('target','var') || isempty(full), target=0; end

if isstruct(A), rp=A.rp; ci=A.ci; 
else [rp ci]=sparse_to_csr(A); 
end

n=length(rp)-1; 
d=-1*ones(n,1); dt=-1*ones(n,1); pred=zeros(1,n);
sq=zeros(n,1); sqt=0; sqh=0; % search queue and search queue tail/head

% start bfs at u
sqt=sqt+1; sq(sqt)=u; 
t=0;  
d(u)=0; dt(u)=t; t=t+1; pred(u)=u;
while sqt-sqh>0
    sqh=sqh+1; v=sq(sqh); % pop v off the head of the queue
    for ri=rp(v):rp(v+1)-1
        w=ci(ri);
        if d(w)<0
            sqt=sqt+1; sq(sqt)=w; 
            d(w)=d(v)+1; dt(w)=t; t=t+1; pred(w)=v; 
            if w==target, return; end
        end
    end
   end
end

function [dCore,D,kdtree] = compute_core_distances( X,k )
    % [dCore,D,kdtree] = compute_core_distances( X,k )
    %
    % computes the core distances of each point in X as:
    %   
    %       dCore(i) = D(i,j); j=k, the kth neareast neighbor of X(i)
    %
    % Inputs:
    %   X - n x m matrix, n = observations, m = variables (dimensions)
    %
    %   k - scalar specifying the number of nearest neighbors 
    %
    % Outputs:
    %   dCore - n x 1 vector of core distances 
    %
    %   D - n x n matrix of pairwise distances if m > 10, or nan otherwise.
    %       In the case of m <= 10, the k-nearest neighbor search is
    %       performed more efficiently using KD trees, thus the n x n
    %       distance matrix is not necessary 
    %
    %   kdtree - if m <= 10 and n > 100, kdtree is a KDTree object of X, 
    %            allowing for fast nearest-neighbor queries in the future.
    %            Else, kdtree = nan;
    %
    % Written by Jordan Sorokin, 10/6/2017

    [n,m] = size( X );
    
    if k == 1
        dCore = zeros( 1,n );
        D = nan;
        return
    end
    
    %if (m >= 10) || (n < 100)
        kdtree = nan;
        if n > 20e3
            [~,dCore,D] = compute_nearest_neighbors( X,X,k-1 ); % slow but memory conservative
        else
            D = compute_pairwise_dist( X );
            dCore = sort( D,1 );
            dCore = dCore(k,:);
        end
    %else
%         kdtree = createns( X,'nsmethod','kdtree' );
%         [neighbors,dCore] = kdtree.knnsearch( X,'k',k );
%         dCore = dCore(:,end);
%         neighbors = neighbors(:,end);
%         D = sparse( double( neighbors ),1:n,dCore,n,n );
    %end
end

function [neighbors,distances,D] = compute_nearest_neighbors( X,Y,k )
    % [neighbors,distances,D] = compute_nearest_neighbors( X,Y,k )
    %
    % comptues the nearest k neighbor distances for each point in X. 
    % X = n x d matrix, with n = observations, d = variables (dimensions)
    % Y = m x d matrix, with m = observations, d = variables (dimensions)
    %
    % neighbors, distances = n x k matrices (default to loop over X)
    %
    % D is a sparse n x m matrix of the pairwise distances between nearest neighbors.
    % Note that it is only reasonable to output D if k << min(n,m). 

    [n,d] = size( X );
    [m,d2] = size( Y );

    assert( d == d2,'# of dimensions of X and Y must match' );
    assert( all( k < [n,m] ),'requested more neighbors than # of points available' );

    % calculate the ranges for looping
    neighbors = zeros( n,k );
    distances = zeros( n,k );
    [start,stop] = find_blockIndex_range( m,n,1e8 );

    for i = 1:numel( start )

        % euclidean distance
        pts = X(start(i):stop(i),:);
        S = compute_pairwise_dist( pts,Y );

        % nearest neighbors
        [S,idx] = sort( S,2,'ascend' ); % sorts smallest - largest columns for each row
        neighbors(start(i):stop(i),:) = idx(:,2:k+1);
        distances(start(i):stop(i),:) = real( sqrt( S(:,2:k+1) ) );
    end

    % compute the sparse distance matrix
    if nargout > 2
        D = sparse( reshape( double( neighbors' ),n*k,1),...
                    reshape( repmat( 1:m,k,1 ),m*k,1 ),...
                    reshape( double( distances ),n*k,1 ),...
                    n,m );
    end

end

function D = compute_pairwise_dist( X,varargin )
    % D = compute_pairwise_dist( X,(Y) )
    %
    % computes a (fast!) pairwise euclidean distance matrix of the inputs.
    % If more than one matrix is provided, will compute pairwise distances
    % between matrices. Else, computes a symmetric pairwise distance of X
    %
    % Both X and the optional matrix Y must have points along the rows, and
    % dimensions of points along columns.

    if nargin == 1
        d = sum( X.^2,2 ); 
        D = real( sqrt( bsxfun( @plus,d,d' ) - 2*(X * X') ) );
    else
        Y = varargin{1};
        if size( Y,2 ) ~= size( X,2 )
            error( 'X and Y must have equal number of dimensions' );
        end

        D = real( sqrt( bsxfun( @plus, sum( X.^2,2 ),...
            bsxfun( @minus, sum( Y.^2,2 )', 2*(X * Y') ) ) ) );
    end
end


function [start,stop] = find_blockIndex_range( n,m,maxSize )
% [start,stop] = find_indexing_range( n,m,(maxSize) )
%
% find the indexing ranges to perform block looping over points "1:m" for 
% nearest-neighbor, etc. This greatly improves performance (not looping over 
% all points) while avoiding massive (> 1 GB) matrices. Default "maxSize",
% which refers to the maximum matrix size allowed, is 100,000,000 entries.
%
% "n" refers to number of elements in first vector, "m" refers to number of
% elements in second. This is more flexible than simply providing one
% number, as neighborhood graphs etc. may not necessarily be square.

if nargin < 3
    maxSize = 1e8; % 1 GB
end

maxPts = ceil( maxSize / n );
remainder = mod( m,maxPts );
start = 1:maxPts:m;
stop = maxPts:maxPts:m-remainder;
if remainder > 0
    stop = [stop,m];
end

end

function children = get_all_subnodes( u,parents )
    % children = get_all_subnodes( u,parents )
    %
    % recursively finds all subnodes of a branching tree 
    % starting from root node "u"
    %
    % Inputs:
    %   u - the root node to search top-down from
    %
    %   parents - a list of parents of length nNode,
    %             where nNode = total # of nodes in the 
    %             tree, and parents(i) indicates the directly-
    %             connected parent to the ith node
    %
    %       i.e. for a branching tree such as:
    %
    %                       (1)
    %                      __|__
    %                    (2)   (3)
    %                   __|__  
    %                 (4)   (5)
    %
    %       the vector 'parent' would be: [0 1 1 2 2]
    %
    % Written by Jordan Sorokin, 10/11/2017
    
    % find children of this node
    subChildren = find(parents == u); 
    if isempty( subChildren )
        children = [];
        return
    end
    
    % add children to growing list
    children = subChildren;
    
    % recursively find children for each new child
    for child = subChildren
        newChildren = get_all_subnodes( child,parents );
        children = [children,newChildren];
    end
end

function [labels,P] = get_cluster_probability( clusters,lambdaMax,coreLambda )
    % [labels,P] = get_cluster_probability( clusters,lambdaMax,coreLambda )
    %
    % computes the probability of point i belonging to cluster j, if 
    % point i is associated with cluster j (if lambdaMax(i,j) > 0)
    %
    % Inputs:
    %   clusters - 1 x K vector of cluster IDs
    %
    %   lambdaMax - n x P >= K matrix of lambda values of point i for
    %               cluster j.  
    %
    %   coreLambda - 1 x K vector of maximum lambda associated with cluster
    %                i or any of its subclusters
    %
    % Outputs:
    %   labels - n x 1 vector of cluster assignments
    %
    %   P - n x 1 vector of assignment probabilities
    %
    % Written by Jordan Sorokin, 10/15/2017

    n = size( lambdaMax,1 );
    P = zeros( n,1 );
    labels = zeros( n,1,'uint8' );
    
    % loop over clusters
    for k = 1:numel( clusters )
        pts = lambdaMax(:,clusters(k)) > 0;
        labels(pts) = k;
        P(pts) = lambdaMax(pts,clusters(k)) / coreLambda(k);
    end
end

function [corePoints,coreLambda] = get_core_points( parents,bestClusters,lambdaMax )
    % [corePoints,coreLambda] = get_core_points( clusters,parents,bestClusters,lambdaMax ) 
    %
    % finds the core points for the final kept clusters by recursively 
    % searching for all subchildren of each kept cluster, and finding which
    % points persisted the longest (maximum lambda value)
    %
    % Inputs: 
    %   parents - 1 x K vector of parents of clusters 1:K
    %
    %   bestClusters - 1 x p <= K vector of kept clusters, as determined by
    %                  the function OFCH.m
    %
    %   lambdaMax - n x K matrix of lambda values for each point i and
    %               cluster j. LambdaMax is the lambda value at which point
    %               i "falls out of" cluster j
    %
    % Outputs:
    %   corePoints - 1 x p cell array of core points, with each cell representing
    %                one cluster ID from the list of best clusters
    %
    %   coreLambda - 1 x p vector of lambda max values associated with the
    %                p best clusters
    %
    % Written by Jordan Sorokin, 10/15/2017
    
    p = length( bestClusters );
    corePoints = cell( 1,p );
    coreLambda = zeros( 1,p );
    
    for k = 1:numel( bestClusters )
        thisClust = bestClusters(k);
        
        % find subnodes of bestCluster(k)
        children = get_all_subnodes( thisClust,parents );
        
        % get maximum lambda value of this cluster or its children
        maxLambda = max( lambdaMax(:,[thisClust,children]),[],2 ); % finds maximum lambda bewteen subnodes
        coreLambda(k) = max( maxLambda ); % finds the maximum lambda across all points and subnodes
        
        % find the core points of this cluster/subclusters
        corePoints{k} = find( maxLambda == coreLambda(k) );
    end
end

function model = hdbscan_fit( X,varargin )
    % model = hdbscan_fit( X,(varargin) )
    %
    % models the rows in X using the heirarchical DBSCAN, as
    % described in the manuscript:
    %
    %       Campello et al. 2015, Hierarchical Density Estimatse for Data
    %       Clustering, Visualization, and Outlier Detection
    %
    % Inputs:
    %   X - an n x m matrix with n = observations, m = variables (dimensions)
    %
    %   (minpts) - scalar specifying the minimum number of points necessary within the
    %            neighborhood radius ball of any point X(i) before that point is considered noise 
    %            for a given eps. radius. (default = 5)
    %
    %   (minclustsize) - minimum number of points for a cluster to be considered a true cluster
    %                    (default = minpts)
    %
    %   (minClustNum) - the minimum # of clusters; the first minClustNum
    %                   parent clusters will have stabilities set to 0
    %                   (default = 1)
    %
    %   (dEps) - positive int used for determining the number of epsilon
    %            (edge weights) to loop over as:
    %
    %               weights = edgeWeights(1:dEps:end)
    %               
    %            default is to use all values of 'edgeWeight'. values of
    %            dEps > 1 result in an increase in model convergence speed
    %            (as a linear scaling) at the expense of rougher
    %            cluster splitting / assignment
    %          
    % Outputs:
    %   model - structure produced by the clustering algorithm containing the following fields:
    %
    %       clusterTree :   a structure with the cluster tree created by the data in X. 
    %                       The tree is condensed in that spurious components (where
    %                       # pts < minclustsize) are eliminated. The
    %                       resulting tree is orders of magnitude smaller
    %                       than the full hierarchy. 
    %
    %                       The structure contains:
    %                           clusters - vector of cluster labels
    %                           parents - their parents
    %                           lambdaMin - minimum lambda (max epsilon) at which each was created
    %                           stability - their stabilities
    %                       
    %                       stabilities are measured as:
    %
    %                           S(j) = 1 - { eps_max(xi) / eps(xi) }
    %
    %       lambda      :   vector of 1 / epsilon radii
    %
    %       lambdaNoise :   the lambda associated with point i beyond which that
    %                       point becomes noise (ID = 0). This is the
    %                       smallest radii (largest lambda) that point i belongs 
    %                       to any cluster j
    %
    %       lambdaMax   :   the lambda associated with point i and cluster
    %                       j after which point i no longer belongs to
    %                       cluster j ... (1 / eps(xi) )
    %
    %       lastClust   :   the final cluster j that point i belonged to
    %                       before becoming noise
    %   
    %       dCore       :   core distances of each X(i) as the distance to
    %                       the 'minpts' nearest neighbor of X(i)
    %
    % Written by Jordan Sorokin
    % 10/5/2017


    %% GLOBALS
    [n,~] = size( X );    
    p = check_inputs( varargin );
    minclustsize = p.minclustsize;     
    minpts = p.minpts;
    dEps = p.dEps;
    maxClustNum = 1000;
    minClustNum = p.minClustNum;
    clear p d

    %% CREATE CONENCTED TREE
    % (a) compute the core distances & mutual reachability for each X(i)    
    [dCore,D] = compute_core_distances( X,minpts );
    D = mutual_reachability( D,dCore );
    
    % (b) create the minimum spanning tree and add self loops
    [nodes(:,1),nodes(:,2),weights] = mst_prim( D ); clear mr
    if ~isa( weights(1),'double' ) || ~isa( weights(1),'logical' )
        weights = double( weights );
        nodes = double( nodes );
    end
    nodes = [nodes;repmat( (1:n)',1,2 )]; weights = [weights;dCore']; % adds self loops
    mst = sparse( [nodes(:,1);nodes(:,2)],[nodes(:,2);nodes(:,1)],[weights;weights],n,n ); % makes it symmetric
    
    % (c) get sorted weight vector for the loop
    epsilon = single( sort( unique( weights ),'descend' ) ); % sorted edge weights (larger = further apart)
    epsilon = epsilon(1:dEps:end);
    nEpsilon = numel( epsilon );

    % (d) pre-allocate our matrices for storage  
    lambdaNoise = zeros( n,1,'single' );            % keeps track of epsilon when X(i) becomes noise
    lastClust = zeros( n,1,'uint32' );              % keepst track of C(j) when X(i) becomes noise
    parentClust = zeros( 1,maxClustNum,'uint32' );  % keeps track of which parent spawned each cluster
    lambdaMin = zeros( 1,maxClustNum,'single' );    % keeps track of max Eps for clust C(j) to appear
    lambdaMax = zeros( n,maxClustNum,'single' );    % keeps track of min Eps for point X(i) still in cluster C(j)
    currentMaxID = uint32( 1 );                     % keeps track of max ID for updating labels of new components
    lambdaMin(1) = 1./epsilon(1);
    newID = ones( n,1,'uint32' );       
    
    %% HIERARCHICAL SEARCH
    for i = 2:nEpsilon
        
        oldID = newID;
        
        % (e) find edges greater than current epsilon value
        idx = weights > epsilon(i);
        if ~any( idx )
            continue
        end
        
        % (f) find the nodes of the cut edges and 
        % remove bad ones (those previously labeled as noise)
        endNodes = nodes(idx,:);
        nodes(idx,:) = [];
        weights(idx) = [];
        for q = 1:nnz( idx )
            mst(endNodes(q,1),endNodes(q,2)) = 0; % cut the edges
            mst(endNodes(q,2),endNodes(q,1)) = 0;
        end
        
        % remove noise
        selfCut = (endNodes(:,1) == endNodes(:,2));
        if any( selfCut )
            newID(endNodes(selfCut,1)) = 0;
            endNodes(selfCut,:) = [];
            if isempty( endNodes )
                continue
            end
        end
        
        % (g) remove noisy nodes and skip loop if no remaining nodes 
        uniqueCuts = unique( endNodes );
        badIDX = (newID(uniqueCuts) == 0);     
        if any( badIDX )
            if all( badIDX )
                continue
            end
            
            % if only some nodes noisy, remove those
            badNodes = uniqueCuts(badIDX);
            endNodes(any( ismember( endNodes',badNodes ) )',:) = [];
        end        

        nCutEdges = size( endNodes,1 );
        
        %% FIND SUB-TREES IN FOREST    
        for k = 1:nCutEdges
            
            % (h) get the connected components from the end nodes
            parent = oldID(endNodes(k,1));
            
            [~,~,temp1] = bfs( mst,endNodes(k,1) );
            [~,~,temp2] = bfs( mst,endNodes(k,2) );
            subTree1 = temp1 > 0;
            subTree2 = temp2 > 0;
            validTree = [nnz( subTree1 ),nnz( subTree2 )] >= minclustsize;            
            
            % (i) check for noise or splits
            % (i.1) - one or both trees are too small
            if ~all( validTree )
                
                if validTree(1)   
                    isNoise = subTree2;
                elseif validTree(2)
                    isNoise = subTree1;
                else
                    isNoise = subTree1 | subTree2;
                end
                
                lastClust(isNoise) = oldID(isNoise);
                lambdaNoise(isNoise) = epsilon(i-1);
                lambdaMax(isNoise,parent) = epsilon(i-1);
                newID(isNoise) = 0;
                
                if ~any( newID )
                    break
                end 
                
            % (i.2) both subtrees valid
            else
                newMax = currentMaxID + 2;
                temp = newMax - 1;
                newID(subTree1) = temp;
                newID(subTree2) = newMax; 
                parentClust(temp:newMax) = parent;
                currentMaxID = newMax;
                lambdaMax(subTree1 | subTree2,parent) = epsilon(i-1);
                lambdaMin([temp,newMax]) = epsilon(i);  
            end
        end
    end
    
    %% COMPUTE CLUSTER STABILITY % OUTPUT MODEL 
    
    % (j) get the condensed tree
    lambdaMin = 1 ./ lambdaMin(1:currentMaxID);
    lambdaMax = 1 ./ lambdaMax(:,1:currentMaxID);
    lambdaMax(isinf(lambdaMax)) = 0;
    parentClust = parentClust(1:currentMaxID);
    
    % (k) compute stabilities
    inClust = ones( n,currentMaxID );
    inClust(lambdaMax == 0) = 0; 
    S = sum( lambdaMax - bsxfun( @times,lambdaMin,inClust ) ); % [EQ 3] ... cluster stability    
    
    % set 1:minClustNum-1 parent clust stabilities to 0
    if minClustNum > 1
        [uniqueParents,parentLocation] = unique( parentClust(parentClust>0) );
        [~,idx] = sort( parentLocation );
        uniqueParents = uniqueParents( idx ); % order in which the parents were split
        nchildren = zeros( 1,numel( uniqueParents ) );
        for i = 1:numel( uniqueParents )
            nchildren(i) = nnz( parentClust == uniqueParents(i) );
        end
        
        % set stabilities for the first stability-ordered N parents that, cumulatively, have
        % minClustNum # of children equal to 0 
        removeParents = uniqueParents( 1:find( cumsum( nchildren ) >= minClustNum,1 ) );
        S(removeParents) = 0;
    end
    
    % Model output
    clusterTree = struct( 'clusters',1:currentMaxID,'parents',parentClust,...
        'lambdaMin',lambdaMin,'stability',S ); 
    
    model = struct( 'lambda',1./epsilon(1:i),'clusterTree',clusterTree,'dCore',dCore,...
                    'lambdaNoise',1./lambdaNoise,'lastClust',lastClust,'lambdaMax',sparse( double( lambdaMax ) ) );

    %% FUNCTIONS
    function p = check_inputs( inputs )
        names = {'minpts','minclustsize','minClustNum','dEps'};
        defaults = {5,nan,1,1};

        p = inputParser;
        for j = 1:numel( names )
            p.addParameter( names{j},defaults{j} );
        end

        parse( p,inputs{:} );
        p = p.Results;

        % change minclustsize if nan
        if isnan( p.minclustsize )
            p.minclustsize = p.minpts;
        end
        
        % check dEps / round to nearest int
        p.dEps = round( p.dEps ) * sign( p.dEps );
        if p.dEps == 0
            p.dEps = 1;
        end
    end
end

function varargout=mst_prim(A,full,u)
% MST_PRIM Compute a minimum spanning tree with Prim's algorithm
% 
% T = mst_prim(A) computes a minimum spanning tree T using Prim's algorithm
% for the spanning tree of a graph with non-negative edge weights.
%
% T = mst_prim(A,0) produces an MST for just the component at A containing
% vertex 1.  T = mst_prim(A,0,u) produces the MST for the component
% containing vertex u.
%
% [ti tj tv] = mst_prim(...) returns the edges from the matrix and does not
% convert to a sparse matrix structure.  This saves a bit of work and is
% required when there are 0 edge weights.
%
% Example:
%   load_gaimc_graph('airports'); % A(i,j) = negative travel time
%   A = -A; % convert to travel time.
%   A = max(A,A'); % make the travel times symmetric
%   T = mst_prim(A);
%   gplot(T,xy); % look at the minimum travel time tree in the US

% David F. Gleich
% Copyright, Stanford University, 2008-2009

%  History:
%  2009-05-02: Added example

% TODO: Add example

if ~exist('full','var') || isempty(full), full=0; end
if ~exist('target','var') || isempty(full), u=1; end

if isstruct(A) 
    rp=A.rp; ci=A.ci; ai=A.ai; 
    check=0;
else
    [rp ci ai]=sparse_to_csr(A); 
    check=1;
end
if (check && any(ai) <= 0), error('gaimc:prim', ...
        'prim''s algorithm cannot handle negative edge weights.'); 
end
if check && ~isequal(A,A'), error('gaimc:prim', ...
        'prim''s algorithm requires an undirected graph.'); 
end
nverts=length(rp)-1; 
d=Inf*ones(nverts,1); T=zeros(nverts,1); L=zeros(nverts,1);
pred=zeros(1,length(rp)-1);

% enter the main dijkstra loop
for iter=1:nverts
    if iter==1, root=u; 
    else 
        root=mod(u+iter-1,nverts)+1; 
        if L(v)>0, continue; end 
    end
    n=1; T(n)=root; L(root)=n; % oops, n is now the size of the heap
    d(root) = 0;
    while n>0
        v=T(1); L(v)=-1; ntop=T(n); T(1)=ntop; n=n-1;
        if n>0, L(ntop)=1;  end         % pop the head off the heap
        k=1; kt=ntop;                   % move element T(1) down the heap
        while true
            i=2*k; 
            if i>n, break; end          % end of heap
            if i==n, it=T(i);           % only one child, so skip
            else                        % pick the smallest child
                lc=T(i); rc=T(i+1); it=lc;
                if d(rc)<d(lc), i=i+1; it=rc; end % right child is smaller
            end
            if d(kt)<d(it), break;     % at correct place, so end
            else T(k)=it; L(it)=k; T(i)=kt; L(kt)=i; k=i; % swap
            end
        end                             % end heap down

        % for each vertex adjacent to v, relax it
        for ei=rp(v):rp(v+1)-1            % ei is the edge index
            w=ci(ei); ew=ai(ei);          % w is the target, ew is the edge weight
            if L(w)<0, continue; end      % make sure we don't visit w twice
            % relax edge (v,w,ew)
            if d(w)>ew
                d(w)=ew; pred(w)=v;
                % check if w is in the heap
                k=L(w); onlyup=0; 
                if k==0
                    % element not in heap, only move the element up the heap
                    n=n+1; T(n)=w; L(w)=n; k=n; kt=w; onlyup=1;
                else kt=T(k);
                end
                % update the heap, move the element down in the heap
                while 1 && ~onlyup,
                    i=2*k; 
                    if i>n, break; end          % end of heap
                    if i==n, it=T(i);           % only one child, so skip
                    else                        % pick the smallest child
                        lc=T(i); rc=T(i+1); it=lc;
                        if d(rc)<d(lc), i=i+1; it=rc; end % right child is smaller
                    end
                    if d(kt)<d(it), break;      % at correct place, so end
                    else T(k)=it; L(it)=k; T(i)=kt; L(kt)=i; k=i; % swap
                    end
                end
                % move the element up the heap
                j=k; tj=T(j);
                while j>1,                       % j==1 => element at top of heap
                    j2=floor(j/2); tj2=T(j2);    % parent element
                    if d(tj2)<d(tj), break;      % parent is smaller, so done
                    else                         % parent is larger, so swap
                        T(j2)=tj; L(tj)=j2; T(j)=tj2; L(tj2)=j; j=j2;
                    end
                end  
            end
        end
    end
    if ~full, break; end
end

nmstedges=0;
for i=1:nverts
    if pred(i)>0, nmstedges=nmstedges+1; end
end
ti = zeros(nmstedges,1); tj=ti; tv = zeros(nmstedges,1);
k=1;
for i=1:nverts
    if pred(i)>0, 
        j = pred(i);
        ti(k)=i; tj(k)=j;
        for rpi=rp(i):rp(i+1)-1
            if ci(rpi)==j, tv(k)=ai(rpi); break; end
        end
        k=k+1;
    end
end
if nargout==1,
    T = sparse(ti,tj,tv,nverts,nverts);
    T = T + T';
    varargout{1} = T;
else
    varargout = {ti, tj, tv};
end

end

function D = mutual_reachability( D,coreDist )
    % D = mutual_reachability( D,coreDist )
    %
    % Finds the mutual reachability between points as:
    %
    %		max{ core(i),core(j),D(i,j) }
    %
    % Inputs:
    %	D - n x n distance matrix
    %	coreDist - 1 x n or n x 1 vector of core distances
    %
    % Outputs:
    %	D - n x n mutual reachability graph, edited in place
    %
    % written by Jordan Sorokin, 10/7/2017

    % find where core(i) or core(j) > D(i,j)
    maxCore = bsxfun( @max,coreDist,coreDist' );
    idx = maxCore > D; 
    D(idx) = maxCore(idx); % ~idx locations left as D

end

function [G,h] = plot_cluster_tree( tree )
    % [G,h] = plot_cluster_tree( tree )
    %
    % plots the condensed cluster tree as a connected 
    % graph, based on the relationship between parent clusters
    % and all sub clusters from the parent
    %
    % Inputs:
    %   tree - a cluster tree created by HDBSCAN.fit() or hdbscan_fit
    %
    % Outputs:
    %   G - the graph object created by linking a parent with all 
    %       sub clusters spawning from the parent
    %
    %   h - the handle to the figure
    %
    % Written by Jordan Sorokin, 10/15/2017

    % check inputs
    clusters = tree.clusters;
    parents = tree.parents;
    weights = tree.lambdaMin;
    if parents(1) == 0
        parents = parents(2:end);
        clusters = clusters(2:end);
        weights = [0,weights(2:end)];
    end
   
    % create a graph object and plot
    G = graph( clusters,parents );
    h = G.plot();
    h.YData = weights;
    h.Parent.YDir = 'reverse';
    h.Parent.YLabel.String = '\lambda';
end
    
function [rp ci ai ncol]=sparse_to_csr(A,varargin)
% SPARSE_TO_CSR Convert a sparse matrix into compressed row storage arrays
% 
% [rp ci ai] = sparse_to_csr(A) returns the row pointer (rp), column index
% (ci) and value index (ai) arrays of a compressed sparse representation of
% the matrix A.
%
% [rp ci ai] = sparse_to_csr(i,j,v,n) returns a csr representation of the
% index sets i,j,v with n rows.
%
% Example:
%   A=sparse(6,6); A(1,1)=5; A(1,5)=2; A(2,3)=-1; A(4,1)=1; A(5,6)=1; 
%   [rp ci ai]=sparse_to_csr(A)
%
% See also CSR_TO_SPARSE, SPARSE

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2008-04-07: Initial version
% 2008-04-24: Added triple array input
% 2009-05-01: Added ncol output
% 2009-05-15: Fixed triplet input

error(nargchk(1, 5, nargin, 'struct'))
retc = nargout>1; reta = nargout>2;

if nargin>1
    if nargin>4, ncol = varargin{4}; end
    nzi = A; nzj = varargin{1};
    if reta && length(varargin) > 2, nzv = varargin{2}; end    
    if nargin<4, n=max(nzi); else n=varargin{3}; end
    nz = length(A);
    if length(nzi) ~= length(nzj), error('gaimc:invalidInput',...
            'length of nzi (%i) not equal to length of nzj (%i)', nz, ...
            length(nzj)); 
    end
    if reta && length(varargin) < 3, error('gaimc:invalidInput',...
            'no value array passed for triplet input, see usage'); 
    end
    if ~isscalar(n), error('gaimc:invalidInput',...
            ['the 4th input to sparse_to_csr with triple input was not ' ...
             'a scalar']); 
    end
    if nargin < 5, ncol = max(nzj); 
    elseif ~isscalar(ncol), error('gaimc:invalidInput',...
            ['the 5th input to sparse_to_csr with triple input was not ' ...
             'a scalar']); 
    end
else
    n = size(A,1); nz = nnz(A); ncol = size(A,2);
    retc = nargout>1; reta = nargout>2;
    if reta,     [nzi nzj nzv] = find(A); 
    else         [nzi nzj] = find(A);
    end
end
if retc, ci = zeros(nz,1); end
if reta, ai = zeros(nz,1); end
rp = zeros(n+1,1);
for i=1:nz
    rp(nzi(i)+1)=rp(nzi(i)+1)+1;
end
rp=cumsum(rp);
if ~retc && ~reta, rp=rp+1; return; end
for i=1:nz
    if reta, ai(rp(nzi(i))+1)=nzv(i); end
    ci(rp(nzi(i))+1)=nzj(i);
    rp(nzi(i))=rp(nzi(i))+1;
end
for i=n:-1:1
    rp(i+1)=rp(i);
end
rp(1)=0;
rp=rp+1;

end


        

