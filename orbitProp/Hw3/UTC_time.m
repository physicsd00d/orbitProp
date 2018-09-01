function UTC = UTC_time(day, month, year, seconds)
% only good for year 2012 at the moment, which is a leap year
% also assuming that the day and seconds are already in UTC
% output is in days

RegYear  = [31 28 31 30 31 30 31 31 30 31 30 31];
LeapYear = [31 29 31 30 31 30 31 31 30 31 30 31];

NumDaysPerMonth = LeapYear;

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

if (year ~= 2012)
    error('This function is currently only good for 2012')
end

UTC = UTC + seconds/(60*60*24);





end