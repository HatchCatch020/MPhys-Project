% C. Monaghan
classdef linTools
    % linTools class
    % 
    % properties:
    %   beamline
    %   ML_algorithm
    %   
    %
    % Methods:
    %   calcRespMatML
    %   calcRespMatC
    %   getBPMorbitdy
    %   getBPMvalues_corr
    %   track_getBPMreadings
    %   setQuadFerrors
    %   setQuadAerrors
    %
    
    properties
       beamline;
       ML_Algorithm     = 'mvn';
    end
    methods
       function respmatML = calcRespMatML(this, numobs_, BPMnoise_, dcorrstrengthScale)
            dcorrstrength = dcorrstrengthScale;    % Set the scale of corrector field strengths
            bpmresn       = BPMnoise_; % Set the resolution of the bpms
            numBPMs = numel(this.beamline.bpmlist);
            %if this.disableBPMs, numBPMs = numel(this.beamline.bpmlist) - numel(this.disabledBPMlist); end

            dbpm   = zeros(numBPMs,  numobs_);
            dcorra = zeros(numel(this.beamline.corrlist), numobs_);

            % Generate input data for ML
            for i = 1:numobs_
                % Generate a random set of changes to corrector strengths...
                dcorr       =  dcorrstrength*randn(numel(this.beamline.corrlist),1);

                % ...and find the change in the closed orbit with these strengths
                bpmorbitdy  = getBPMorbitdy(this, dcorr) + randn(numBPMs,1)*bpmresn;

        %         for n = 1:numel(this.beamline.corrlist)
        %            dcorr(n) = dcorr(n)*(this.beamline.corrlist{n}.length/this.beamline.beam.rigidity); 
        %         end

                % Record the observation
                dbpm(:,i)   = bpmorbitdy;
                dcorra(:,i) = dcorr;

            end

            xmat  = dcorra';
            ymat  = dbpm';

            xcell = cell(1,numobs_);

            for i = 1:numobs_
                xcell{i} = [kron([xmat(i,:)],eye(numBPMs))];
            end

            % Fit a response matrix to the observed changes in trajectory resulting
            % from given changes in corrector strengths
            [beta,sigma,E,V] = mvregress(xcell,ymat, 'algorithm', this.ML_Algorithm);

            % Calculate the error on the fit
            se = sqrt(diag(V));

            % beta is the response matrix found from the machine learning approach
            respmatML = reshape(beta, [numBPMs, numel(this.beamline.corrlist)]);
       end % function CalcRespMatML

        function respmatC = calcRespMatC(this, dcorrfieldScale)
            % Make sure all kickers are set to zero
            for i = 1:numel(this.beamline.corrlist)
                this.beamline.corrlist{i}.field = [0, 0];
            end

            % Get the bpm readings with no correctors applied
            bpmdy_nocorr = track_getBPMreadings(this);

            % Prepare the matrix for the response matrix
            respmatC = zeros(numel(this.beamline.bpmlist), numel(this.beamline.corrlist));

            for i = 1:numel(this.beamline.corrlist)
                this.beamline.corrlist{i}.field = [dcorrfieldScale, 0];

                % Get the bpm readings with the new correctors applied
                bpmdy = track_getBPMreadings(this) - bpmdy_nocorr;

                for j = 1:numel(this.beamline.bpmlist)
                    % Calculate the components of the response matrix
                    respmatC(j,i) = (bpmdy(j)/dcorrfieldScale);%*(this.beamline.beam.rigidity/this.beamline.corrlist{i}.length);
                end

                % Turn off this corrector
                this.beamline.corrlist{i}.field = [0, 0];
            end
        end % function CalcRespMatC

        function bpmorbitdy_ = getBPMorbitdy(this, dcorrStrength_)
            % Track the particle and retreive the BPM readings 
            bpmorbitdy = track_getBPMreadings(this);

            % Apply changes in correctors
            for i = 1:numel(this.beamline.corrlist)
               this.beamline.corrlist{i}.field = this.beamline.corrlist{i}.field + [dcorrStrength_(i), 0];
            end
            
            % Get the new BPM readings and calculate the change
            bpmorbitdy = track_getBPMreadings(this) - bpmorbitdy;

            % Return corrector strengths to previous value
            for i = 1:numel(this.beamline.corrlist)
               this.beamline.corrlist{i}.field = this.beamline.corrlist{i}.field - [dcorrStrength_(i), 0];
            end
            bpmorbitdy_ = bpmorbitdy;
        end% function getBPMorbitdy
        
        function bpmreadings_ = getBPMvalues_corr(this, dcorr, bpmresn)
            % Apply changes in correctors
            for i = 1:numel(this.beamline.corrlist)
               this.beamline.corrlist{i}.field = this.beamline.corrlist{i}.field + [dcorr(i), 0];
            end
            
            % Get the new BPM readings
            bpmreadings_ = track_getBPMreadings(this) + randn(numel(this.beamline.bpmlist),1)*bpmresn;

            % Return corrector strengths to previous value
            for i = 1:numel(this.beamline.corrlist)
               this.beamline.corrlist{i}.field = this.beamline.corrlist{i}.field - [dcorr(i), 0];
            end
            
        end % function getBPMvalues_corr

        function bpmreadings_ = track_getBPMreadings(this)
            this.beamline.beam.momentum = this.beamline.momentum0;
            this.beamline.beam.particles = [0 0 0 0 0 0]';
            
            bpmorbitdy = zeros(numel(this.beamline.bpmlist), 1);

            % Reset the BPM buffers
            for j = 1:numel(this.beamline.bpmlist)
                this.beamline.bpmlist{j}.ResetBuffer(1);
            end

            n = length(this.beamline.bl.componentlist);
            this.beamline.bl.Track([1 n],this.beamline.beam);

            % Retreive the BPM readings
            for j = 1:numel(this.beamline.bpmlist)
                bpmorbitdy(j) = this.beamline.bpmlist{j}.buffer(2);    
            end
           
            if this.disableBPMs
                for i=1:numel(this.disabledBPMlist)
                    bpmorbitdy(this.disabledBPMlist(i)) = 0;
                end
            end

            bpmreadings_ = bpmorbitdy;
        end % function track_getBPMreadings

        function setQuadFerrors(this, bool, quadError_)
            if bool
                for n = 1:numel(this.beamline.quadlist)
                    gradient = this.beamline.quadlist{n}.gradient;
                    x = (1 + quadError_) + quadError_*randn;
                    this.beamline.quadlist{n}.gradient = gradient * x;
                end
            elseif bool == false
                for n = 1:numel(this.beamline.quadlist)
                    this.beamline.quadlist{n}.gradient = this.beamline.gradient0(n);
                end
            end
            
        end % function setQuadFerrors
        
        function setQuadAerrors(this, bool, deltaY)
            Lkick       = this.beamline.quaderrlist{1}.length;
            deltaYvals  = randn(numel(this.beamline.quadlist),1)*deltaY;
            if bool 
                for n=1:numel(this.beamline.quadlist)
                   intgField = deltaYvals(n)*this.beamline.quadlist{n}.gradient ...
                                       *this.beamline.quadlist{n}.length;
                   this.beamline.quaderrlist{n}.field = [intgField/Lkick, 0];
                end
            elseif bool == false
                for n = 1:numel(this.beamline.quaderrlist)
                    this.beamline.quaderrlist{n}.field = [0,0];
                end
            end
        end % function setQuadAerrors

    end % methods
    
end % classdef