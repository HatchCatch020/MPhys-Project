% Define the CLARA beamline based on the input of a csv file containing the
% information for each component.

clear

%% Open File
filename = 'CLARA_lattice.csv';
fileID = fopen(filename,'r','n','UTF-8');

CLARA_lattice = textscan(fileID, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]', 'Delimiter', ',', 'TextType', 'string',  'ReturnOnError', false);

fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(CLARA_lattice{1}),length(CLARA_lattice)-1);
for col=1:length(CLARA_lattice)-1
    raw(1:length(CLARA_lattice{col}),col) = mat2cell(CLARA_lattice{col}, ones(length(CLARA_lattice{col}), 1));
end
numericData = NaN(size(CLARA_lattice{1},1),size(CLARA_lattice,2));

for col=[1,4,5,6,7,8,9,10,11,12,13,14,15,16]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = CLARA_lattice{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,4,5,6,7,8,9,10,11,12,13,14,15,16]);
rawStringColumns = string(raw(:, [2,3]));
NumericInfo = cell2mat(rawNumericColumns);

%% Count for each component
% counters = [nDrift nQuadf nQuadd nKicker nMonitor nLcavity nSbend]
c = zeros(1, 7);

for i=1:length(rawStringColumns)
    switch rawStringColumns(i,2)
        case "Drift"
            c(1) = c(1) + 1;
        case "Quadrupole"
            if NumericInfo(i,11) > 0
                c(2) = c(2) + 1;
            else
                c(3) = c(3) + 1;
            end   
        case "Kicker"
            c(4) = c(4) + 1;
        case "Monitor"
            c(5) = c(5) + 1;
%         case "Lcavity" % Temp set to drift for now
%             c(6) = c(6) + 1;
        case "Sbend"
            c(7) = c(7) + 1;
        case {"Ecollimator","Rcollimator","Lcavity"} % set as drift
            c(1) = c(1) + 1;
        otherwise
            continue
    end
end

clearvars -except c NumericInfo rawStringColumns


%% DefineBeamline
beam  = Beam(Electron);
momentum0 = 35.0 * PhysicalUnits.MeV / PhysicalConstants.SpeedOfLight;

beam.momentum = momentum0;

driftlist = cell(1, c(1));
quadflist = cell(1, c(2));
quaddlist = cell(1, c(3));
corrlist  = cell(1, c(4));
bpmlist   = cell(1, c(5));
Lcavitylist = cell(1,c(6));
Sbendlist   = cell(1,c(7));

% Lcavityfreq  = cell(1,c(7));
% Lcavityphase = cell(1,c(7));

pvals    = zeros(length(NumericInfo)+1,1);
pvals(1) = beam.momentum;

bl = Beamline;
MasterOscillator.SetFrequency(2998500000);

c = ones(1, 7);

% For each component in the file execute the corresponding block of code
for n=1:length(NumericInfo)
    switch rawStringColumns(n,2)
        case "Drift"
            driftlist{c(1)}      = Drift;
            driftlist{c(1)}.length = NumericInfo(n,3);
            bl.AppendComponent(driftlist{c(1)});
            c(1) = c(1) + 1;
        case "Quadrupole"
            if NumericInfo(n,11) > 0
               quadflist{c(2)}         = Quadrupole;
               quadflist{c(2)}.length    = NumericInfo(n,3);
               quadflist{c(2)}.gradient  = NumericInfo(n,11) * beam.rigidity;
               bl.AppendComponent(quadflist{c(2)});
               c(2) = c(2) + 1;
            else
               quaddlist{c(3)}         = Quadrupole;
               quaddlist{c(3)}.length    = NumericInfo(n,3);
               quaddlist{c(3)}.gradient  = NumericInfo(n,11) * beam.rigidity;
               bl.AppendComponent(quaddlist{c(3)});
               c(3) = c(3) + 1;
            end
        case "Kicker"
            corrlist{c(4)}      = OrbitCorrector;
            corrlist{c(4)}.length = NumericInfo(n,3);
            corrlist{c(4)}.field  = [0, 0];
            bl.AppendComponent(corrlist{c(4)});
            c(4) = c(4) + 1;
        case "Monitor"
            bpmlist{c(5)} = BeamPositionMonitor;
            bl.AppendComponent(bpmlist{c(5)});
            c(5) = c(5) + 1;
            % Drift added by AW 12/2/2021
            driftlist{c(1)} = Drift;
            driftlist{c(1)}.length = NumericInfo(n,3);
            bl.AppendComponent(driftlist{c(1)});
            c(1) = c(1) + 1;
        case "Lcavity"
            % case "Lcavity" modified by AW 12/2/2021
            Lcavitylist{c(6)} = LinacStructure;
            Lcavitylist{c(6)}.structuretype = 'TravellingWave';
            Lcavitylist{c(6)}.length = NumericInfo(n,3);
            Lcavitylist{c(6)}.ncell = round(3*NumericInfo(n,3)*NumericInfo(n,14)/PhysicalConstants.SpeedOfLight);
            Lcavitylist{c(6)}.harmonic = round(NumericInfo(n,14)/MasterOscillator.GetFrequency());
            Lcavitylist{c(6)}.voltage = -NumericInfo(n,12) * NumericInfo(n,3);
            Lcavitylist{c(6)}.phase = NumericInfo(n,13)*2*pi/360;
            Lcavitylist{c(6)}.globalclock = false;
            bl.AppendComponent(Lcavitylist{c(6)});
            beam.energy = beam.energy + beam.species.charge * Lcavitylist{c(6)}.voltage * cos(Lcavitylist{c(6)}.phase);
            c(6) = c(6) + 1;
        case "Sbend"
            Sbendlist{c(7)} = Dipole;
            Sbendlist{c(7)}.length     = NumericInfo(n,3);
            Sbendlist{c(7)}.curvature  = 1/NumericInfo(n,10);
            theta = NumericInfo(n,3)/NumericInfo(n,10);
            Sbendlist{c(7)}.e1 = theta/2;
            Sbendlist{c(7)}.e2 = theta/2;
            Sbendlist{c(7)}.field      = beam.rigidity * Sbendlist{c(7)}.curvature;
            bl.AppendComponent(Sbendlist{c(7)});
            c(7) = c(7) + 1;
        case {"Ecollimator","Rcollimator","Lcavity"}
            driftlist{c(1)}      = Drift;
            driftlist{c(1)}.length = NumericInfo(n,3);
            bl.AppendComponent(driftlist{c(1)});
            c(1) = c(1) + 1;
        otherwise
            continue
    end
    pvals(numel(bl.componentlist)+1) = beam.momentum;
end

fprintf('Initial beam momentum = %6.5g MeV/c\n', momentum0 / PhysicalUnits.MeV * PhysicalConstants.SpeedOfLight);
fprintf('  Final beam momentum = %6.5g MeV/c\n', beam.momentum / PhysicalUnits.MeV * PhysicalConstants.SpeedOfLight);

beam.momentum = momentum0;

svals = bl.ComputePositions();

fprintf('Total beamline length = %6.5g m\n', svals(end) );

% Setup structure for all useful variables in this script
beamline.driftlist = driftlist;
beamline.quadflist = quadflist;
beamline.quaddlist = quaddlist;
beamline.corrlist  = corrlist;
beamline.bpmlist   = bpmlist;
beamline.Lcavitylist = Lcavitylist;
beamline.Sbendlist   = Sbendlist;
beamline.beam = beam;
beamline.bl = bl;
beamline.momentum0 = momentum0;

% clearvars -except beamline driftlist quadflist quadflist corrlist bpmlist Lcavitylist Sbendlist beam bl
