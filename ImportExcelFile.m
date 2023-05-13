function [TrackN, SliceN, X1, Y1, Distance, Velocity, PixelValue] = ImportExcelFile(FileLocation, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  [TRACKN, SLICEN, X1, Y1, DISTANCE, VELOCITY, PIXELVALUE] =
%  IMPORTFILE(FILE) reads data from the first worksheet in the Microsoft
%  Excel spreadsheet file named FILE.  Returns the data as column
%  vectors.
%
%  [TRACKN, SLICEN, X1, Y1, DISTANCE, VELOCITY, PIXELVALUE] =
%  IMPORTFILE(FILE, SHEET) reads from the specified worksheet.
%
%  [TRACKN, SLICEN, X1, Y1, DISTANCE, VELOCITY, PIXELVALUE] =
%  IMPORTFILE(FILE, SHEET, DATALINES) reads from the specified worksheet
%  for the specified row interval(s). Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  [TrackN, SliceN, X1, Y1, Distance, Velocity, PixelValue] = ImportExcelFile("E:\UC_Merced\2022_Spring\BioPhys\FinalProject\ExperimentalData.xlsx", "Results from MAX_20210920_Sox17", [1, 3432]);
%
%  See also READTABLE.


%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [1, 3432];
end


%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":G" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["TrackN", "SliceN", "X1", "Y1", "Distance", "Velocity", "PixelValue"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% Import the data
tbl = readtable(FileLocation, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":G" + dataLines(idx, 2);
    tb = readtable(FileLocation, opts, "UseExcel", false);
    tbl = [tbl; tb]; %#ok<AGROW>
end


%% Convert to output type
TrackN = tbl.TrackN;
SliceN = tbl.SliceN;
X1 = tbl.X1;
Y1 = tbl.Y1;
Distance = tbl.Distance;
Velocity = tbl.Velocity;
PixelValue = tbl.PixelValue;
end