function seeg_secs = convert_transitions_to_ieeg_secs(st_dates,st_times)

ref_start_time = datetime('01/01/2000 00:00:00','InputFormat','MM/dd/yyyy hh:mm:ss');

new_st_dates = st_dates;
for i = 1:length(new_st_dates) 
    if strcmp(new_st_dates{i},' ') % replace empty with last date
        new_st_dates{i} = new_st_dates{i-1};
    end
end
st_dates = new_st_dates;

% Combine dates and times
dts = cellfun(@(x,y) [x, ' ',y],st_dates,st_times,'uniformoutput',false);

% convert to seconds into ieeg file
dt = cellfun(@(x) datetime(x,'InputFormat','dd-MMM-yyyy HH:mm:ss'),dts);
seeg_secs = seconds(dt-ref_start_time);    

end