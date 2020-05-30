function [ mi_inf ] = mutual_information( data, bins, subpop_fractions )
% MUTUAL_INFORMATION calculates the mutual information between the
% responses (table values) and signals (individual tables) of a set of
% data. This is done according to the paper "Information Transuction
% Capacity of Noisy Biochemical Signaling Networks" by Cheong et al. 2011.
%
% The mutual information is calculated by the formula:
%   I(R;S) = SUM_S( SUM_R( P(R,S)log2( P(R,S)/(P(R)P(S)))
% 
% Where I(R;S) is the mutual information between a Response (R) and a
% Signal (S), P(R,S) is the joint probability of a given R,S value, P(R)
% and P(S) are the marginal probabilities, and SUM_R is the sum over all R
% values.
%
% This value must be calculated by binning both the Signal and the Response
% data, but the choice of binning can have a dramatic effect on the measure
% of I(R;S). This function will vary the Response binning and thus give an
% estimate of unbiased mutual information, but the user will need to vary
% the Signal binning by calling this function on differentially-sampled
% Signal values.
%
%
% Input     -  DATA:            A cell array of vectors, each
%                               corresponding to an input value (S in the
%                               Cheong2011 paper). The contents of the
%                               vectors make up the response events (R) for
%                               a particular input. Lengths can vary.
%                               Multiple S values are allowed per event, so
%                               the dimensionality of each array must be
%                               the same and must be MxN, where each m
%                               corresponds to an event (datapoint) and
%                               each n corresponds to a dimesion (measure).
%           -  BINS:            a vector containing the number of bins to
%                               break the response values (R) into. This is
%                               needed, along with the subpop_fractions, in
%                               order to estimate the unbiased mutual
%                               information. E.g. [10,20,30,40]. Must be
%                               integers > 2.
%           -  SUBPOP_FRACTIONS:a vector containing the fractional size of
%                               each jacknife sample to take when
%                               estimating mutual information (e.g.
%                               [0.5,0.8,1]. Must be in [0,1].
%           
% Output    -  MI_INF:          all I_infinity values gotten by fitting a
%                               regression line to the MI values within
%                               each bin number calculated by varying the
%                               population fraction size. The mean/median
%                               of a subset of these values and the sd/mad
%                               should be taken as the unbiased measure of
%                               mutual information (see SOM Section 2.1 of
%                               the Cheong2011 reference).
%
% Copyright     : 2012 Altschuler/Wu Lab
% Created       : 2012.08.28
% Author        : A.Coster
%
% CHANGELOG
% 20120906 - Incomplete modifications to generalize to N dimensions

    % Check the dimensionality of the data.
    dimensions = size( data{1} ) ;
    dimensions = dimensions(2) ;

    % For each combination of BIN and SUBPOP_FRACTION, a value of mutual
    % information (MI) is calculated. These are all biased measures.
    
    % So we need a matrix to collect the biased MI measures, indexed by
    % (BIN, SUBPOP_FRACTION)
    MI_biased    =  zeros( length(bins), length(subpop_fractions) ) ;
    subpop_sizes =  zeros(1, length(subpop_fractions) ) ;
    
    % Loop through each SUBPOP/BIN combo to populate MI_BIASED
    for subpop = 1:length(subpop_fractions),
        

        for bin = 1:length(bins),
            
            
            % Collect all datapoints to set bins. This must be done by
            % dimension. The P(R1,R2) 
            for dim = 1:length(dimensions)
                all_R   = [] ;
                subdata = data ;

                for S = 1:length(data),
                    d          = subdata{S} ;

                    % Randomly choose the subdata locations (Jacknife)
                    positions  = zeros(1,length(d)) ;
                    positions(1:round(subpop_fractions(subpop)*length(d))) = true ;
                    positions = randsample( positions, length(positions), false ) ;

                    subdata{S} = d(logical(positions)) ;

                    all_R = [all_R d(logical(positions)) ] ;
                end

                subpop_sizes( subpop ) = length(all_R) ;

                % Now we have all of the data, and can go through the process of
                % re-binning the Response values according to BINS.

                % Bin the data into 1/bin bins with roughly equal numbers of
                % events in each.
                bin_values = quantile( all_R, (0:bins(bin))/bins(bin)) ;
                bin_values(end) = bin_values(end) + 1 ;

                % Get the counts for each bin. These will end up being the P(R)
                % values after dividing by the total.
                PR = histc( all_R, bin_values ) ;
                PR = PR(1:end-1) ;
                PR = PR/sum(PR) ;

                PS = cellfun( @length, subdata ) ;
                PS = PS/sum(PS) ;

                % All that remains are the P(R,S) values
                PRS = zeros( bins(bin), length(subdata) ) ;

                for S = 1:length(subdata),
                    d        = subdata{S} ;
                    h        = histc( d, bin_values ) ;
                    PRS(:,S) = h(1:end-1) ;
                end

                PRS = PRS/sum(PRS(:)) ;
            end
            
            % Everything for calculating MI for this bin/subdata pair is
            % complete.
            MI  = PRS .* log2(PRS ./ (repmat(PS,length(PR),1).*repmat(PR',1,length(PS)))) ;
            MI  = MI(~isnan(MI)) ;
            MI_biased( bin, subpop ) = sum(MI(:)) ;
            
        end
    end
    
    % The biased measures of MI are now collected, and so it is time to
    % estimate the unbiased measure. This is done by fitting a line to the
    % points within each bin number (each point being from a specific
    % sample size).
    mi_inf = zeros(1,length(bins)) ;
    for bin = 1:length(bins)
        
        y = MI_biased( bin, : ) ;
        x = 1 ./ subpop_sizes ;
        p = polyfit( x, y, 1 )  ;
        mi_inf( bin ) = p(2) ;
    end
end

