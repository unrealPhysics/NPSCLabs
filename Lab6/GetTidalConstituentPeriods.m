function PeriodDays = GetTidalConstituentPeriods(NConstit)
%Returns a vector contrining the periods of the specified number of tidal constituents
%up to the maximum of 37
%Data is from the NOAA web site
%
%Alec Duncan, Curtin University, 6/4/2018

ConstituentPeriodHours = [12.42060122 12.00000000 12.65834824 23.93446966 6.21030061 ...
    25.81934166 4.14020040 8.17713995 6.00000000 6.26917391 12.62600441 4.00000000 ...
    12.87175763 12.90537448 22.30607420 12.22177415 24.00000000 24.83324836 23.09847677 ...
    661.30920485 4382.90520872 8765.82108959 354.36705221 327.85896891 26.72305330 ...
    26.86835667 12.01644920 11.98359578 28.00622255 24.06589016 11.60695156 8.28040081 ...
    12.19162021 8.38630297 11.96723479 3.10515031 6.10333928];

NAvail = length(ConstituentPeriodHours);
if NConstit > NAvail
    NConstit = NAvail;
end

PeriodDays = ConstituentPeriodHours(1:NConstit)/24;