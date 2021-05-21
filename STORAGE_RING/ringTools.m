% C. Monaghan
classdef ringTools
    % ringTools class
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
        end % function calcRespMatML
        
        function respmatC=calcRespMatC(this, dcorrfieldScale)
            for j = 1:numel(this.beamline.bpmlist)
                this.beamline.bpmlist{j}.ResetBuffer(1);
            end

            for i = 1:numel(this.beamline.corrlist)
                this.beamline.corrlist{i}.field = [0, 0];
            end

            % compute the closed orbit without changes in corr strength
            closedorbit_nocorr = ComputeClosedOrbit(this.beamline.bl, this.beamline.beam);

            this.beamline.beam.particles = closedorbit_nocorr(:,1);
            this.beamline.bl.Track([1, numel(this.beamline.bl.componentlist)], this.beamline.beam);

            % Prepare the array BPM values without correctors
            bpmorbity_nocorr = zeros(numel(this.beamline.bpmlist), 1);

            % fill the array with the values
            for j=1:numel(this.beamline.bpmlist)
               bpmorbity_nocorr(j) =  this.beamline.bpmlist{j}.buffer(2);
            end

            % Prepare the matrix for the response matrix
            respmatC = zeros(numel(this.beamline.bpmlist), numel(this.beamline.corrlist));

            % Main loop for applying the changes in corrector strength. Turning on
            % a different corr in each iteration.
            for i=1:numel(this.beamline.corrlist)

                % Set the corr strength as the input of the function.
                dcorrfield = dcorrfieldScale;

                this.beamline.corrlist{i}.field = [dcorrfield, 0];

                % Compute the closed orbit with this change in corrector strength 
                closedorbit = ComputeClosedOrbit( this.beamline.bl,  this.beamline.beam);

                for j = 1:numel(this.beamline.bpmlist)
                     this.beamline.bpmlist{j}.ResetBuffer(1);
                end

                 this.beamline.beam.particles = closedorbit(:,1);
                 this.beamline.bl.Track([1, numel( this.beamline.bl.componentlist)], this.beamline.beam);

                for j=1:numel( this.beamline.bpmlist)
                    % Calculate the change in BPM as a result of this change in
                    % corr strength.
                    dy =  this.beamline.bpmlist{j}.buffer(2) -  bpmorbity_nocorr(j);
                    % Calculate this component of the response matrix 
                    respmatC(j,i) = dy/dcorrfield;
                end

                % Turn off this corrector in preperation for the next itteration
                this.beamline.corrlist{i}.field = [0, 0];

            end

        end % function calcRespMatC
        
        function bpmorbitdy_ = getBPMorbitdy(this, dcorrStrength_)

            bpmorbitdy = zeros(numel(this.beamline.bpmlist),1);

            % Find the closed orbit at the start
            closedorbit = ComputeClosedOrbit(this.beamline.bl, this.beamline.beam);

            for j = 1:numel(this.beamline.bpmlist)
                this.beamline.bpmlist{j}.ResetBuffer(1);
            end

            this.beamline.beam.particles = closedorbit(:,1);
            this.beamline.bl.Track([1, numel(this.beamline.bl.componentlist)], this.beamline.beam);

            for j = 1:numel(this.beamline.bpmlist)
                bpmorbitdy(j) = this.beamline.bpmlist{j}.buffer(2);
            end

            % Now apply the changes in corrector strength
            for i = 1:numel(this.beamline.corrlist)
               this.beamline.corrlist{i}.field = this.beamline.corrlist{i}.field + [dcorrStrength_(i), 0];
            end

            % Find the change in closed orbit resulting from the change in
            % corrector strengths
            closedorbit = ComputeClosedOrbit(this.beamline.bl,this.beamline.beam);

            for j = 1:numel(this.beamline.bpmlist)
                this.beamline.bpmlist{j}.ResetBuffer(1);
            end

            this.beamline.beam.particles = closedorbit(:,1);
            this.beamline.bl.Track([1, numel(this.beamline.bl.componentlist)], this.beamline.beam);

            for j = 1:numel(this.beamline.bpmlist)
                bpmorbitdy(j) = this.beamline.bpmlist{j}.buffer(2) - bpmorbitdy(j);
            end

            % Put the corrector strengths back where they were
            for i = 1:numel(this.beamline.corrlist)
               this.beamline.corrlist{i}.field = this.beamline.corrlist{i}.field - [dcorrStrength_(i), 0];
            end

            % Finally, return the result
            bpmorbitdy_ = bpmorbitdy;

        end % function getBPMorbitdy
        
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
            bpmorbit = zeros(numel(this.beamline.bpmlist),1);

            % Find the closed orbit at the start
            closedorbit = ComputeClosedOrbit(this.beamline.bl, this.beamline.beam);

            for j = 1:numel(this.beamline.bpmlist)
                this.beamline.bpmlist{j}.ResetBuffer(1);
            end

            this.beamline.beam.particles = closedorbit(:,1);
            this.beamline.bl.Track([1, numel(this.beamline.bl.componentlist)], this.beamline.beam);

            for j = 1:numel(this.beamline.bpmlist)
                bpmorbit(j) = this.beamline.bpmlist{j}.buffer(2);
            end

            bpmreadings_ = bpmorbit;
        end % function track_getBPMreadings
        
        function setQuadFerrors(this, bool, quadError_)
            if bool
                for n = 1:numel(this.beamline.quadlist)
                    gradient = this.beamline.quadlist{n}.gradient;
                    this.beamline.quadlist{n}.gradient = gradient * (1 + quadError_*randn);
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
end % classdef ringTools