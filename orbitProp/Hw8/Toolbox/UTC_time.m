function UTC = UTC_time(day, month, year, hours, minutes, seconds, offset)
% Updated May 28 2012
% Good for any year
% offset is in hours
% output is in days

RegYear  = [31 28 31 30 31 30 31 31 30 31 30 31];
LeapYear = [31 29 31 30 31 30 31 31 30 31 30 31];

if (mod(abs(year-1600),4) == 0)
    NumDaysPerMonth = LeapYear;
else
    NumDaysPerMonth = RegYear;
end

UTC = 0;
i = 1;
while (i < month)
    UTC = UTC + NumDaysPerMonth(i);
    i = i+1;
end

if (day <= NumDaysPerMonth(month))
    UTC = UTC + day;
else
    error('Not that many days in the month')
end

% if (year ~= 2012)
%     error('This function is currently only good for 2012')
% end

UTC = UTC + (hours+offset)/24 + minutes/(24*60) + seconds/(60*60*24);





end