function [DayNum, TideHeight] = ReadQueenslandTideData(FileName, RefYear)
%This function reads "csv" historical tide data downloaded from the Queensland
%government website:
%https://www.msq.qld.gov.au/Tides/Open-data
%
%Usage: 
%[DayNum, TideHeight] = ReadQueenslandTideData(FileName, RefYear);
%FileName is a string containing the name of the file to be read (including path and extension)
%RefYear - A DayNum of 0 corresponds to 00:00:00 on the first of January of
%this year (numeric)
%
%DayNum is a column vector of Matlab date numbers (see help on datenum function)
%TideHeight is a corresponding column vector of tide heights

AllDat = load(FileName, '-ascii');

TideHeight = AllDat(:, 2);

%Some fancy footwork is required to convert the (numeric) date-time into separate variables
%and then into Matlab datenums
DateTime = AllDat(:, 1);

Minute = rem(DateTime, 100);
DateTime = floor(DateTime/100);

Hour = rem(DateTime, 100);
DateTime = floor(DateTime/100);

Year = rem(DateTime, 10000);
DateTime = floor(DateTime/10000);

Month = rem(DateTime, 100);
Day = floor(DateTime/100);

NPt = length(Minute);

%Build into an array of date vectors
AllDateVec = [Year, Month, Day, Hour, Minute, zeros(NPt, 1)];

DayNum = datenum(AllDateVec)-datenum([RefYear, 1, 1, 0, 0, 0]);
