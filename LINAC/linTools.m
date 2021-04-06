classdef linTools
    properties
       beamline;
       ML_Algorithm = 'mvn';
    end
    methods
       function respmatML = calcRespMatML(this, numobs_, BPMnoise_, dcorrstrengthScale)
            dcorrstrength = dcorrstrengthScale;    % Set the scale of corrector field strengths
            bpmresn       = BPMnoise_; % Set the resolution of the bpms

            dbpm   = zeros(numel(this.beamline.bpmlist),  numobs_);
            dcorra = zeros(numel(this.beamline.corrlist), numobs_);

            % Generate input data for ML
            for i = 1:numobs_

                % Generate a random set of changes to corrector strengths...
                dcorr       =  dcorrstrength*randn(numel(this.beamline.corrlist),1);

                % ...and find the change in the closed orbit with these strengths
                bpmorbitdy  = getBPMorbitdy(this, dcorr) + randn(numel(this.beamline.bpmlist),1)*bpmresn;

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
                xcell{i} = [kron([xmat(i,:)],eye(numel(this.beamline.bpmlist)))];
            end

            % Fit a response matrix to the observed changes in trajectory resulting
            % from given changes in corrector strengths
            [beta,sigma,E,V] = mvregress(xcell,ymat, 'algorithm', this.ML_Algorithm);

            % Calculate the error on the fit
            se = sqrt(diag(V));

            % beta is the response matrix found from the machine learning approach
            respmatML = reshape(beta, [numel(this.beamline.bpmlist), numel(this.beamline.corrlist)]);
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
            bpmorbitdy_ = track_getBPMreadings(this) - bpmorbitdy;

            % Return corrector strengths to previous value
            for i = 1:numel(this.beamline.corrlist)
               this.beamline.corrlist{i}.field = this.beamline.corrlist{i}.field - [dcorrStrength_(i), 0];
            end

        end% function getBPMorbitdy
        
        function bpmreadings_ = getBPMvalues_corr(this, dcorr, bpmresn)
            % Apply changes in correctors
            for i = 1:numel(this.beamline.corrlist)
               this.beamline.corrlist{i}.field = this.beamline.corrlist{i}.field + [dcorr(i), 0];
            end
            
            % Get the new BPM readings and calculate the change
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

            % Track the particle through the beamline
            for n = 1:length(this.beamline.bl.componentlist)
                this.beamline.bl.Track([n n],this.beamline.beam);
            end

            % Retreive the BPM readings
            for j = 1:numel(this.beamline.bpmlist)
                bpmorbitdy(j) = this.beamline.bpmlist{j}.buffer(2);    
            end

            bpmreadings_ = bpmorbitdy;
        end % function track_getBPMreadings

        function setQuadFerrors(this, bool, quadError_)
            if bool
                for n = 1:numel(this.beamline.quadflist)
                    fgradient = this.beamline.quadflist{n}.gradient;
                    this.beamline.quadflist{n}.gradient = fgradient * (1 + quadError_*randn);
                end
                for n = 1:numel(this.beamline.quaddlist)
                    dgradient = this.beamline.quaddlist{n}.gradient;
                    this.beamline.quaddlist{n}.gradient = dgradient * (1 + quadError_*randn);
                end
            elseif bool == false
                for n = 1:numel(this.beamline.quadflist)
                    this.beamline.quadflist{n}.gradient = this.beamline.fgradient0(n);
                end
                for n = 1:numel(this.beamline.quaddlist)
                    this.beamline.quaddlist{n}.gradient = this.beamline.dgradient0(n);
                end
            end
            
        end % function setQuadFerrors
        
        function setQuadAerrors(this, bool, quaderrorstrength)
            if bool 
                for n = 1:numel(this.beamline.quaderrlist)
                    this.beamline.quaderrlist{n}.field = [quaderrorstrength(n),0];
                end
            else  
                for n = 1:numel(this.beamline.quaderrlist)
                    this.beamline.quaderrlist{n}.field = [0,0];
                end
            end
        end % function setQuadAerrors
        
    end % methods
    
end % classdef