function pt = find_sleep_wake_annotations(pt)

for p = 1:length(pt)
    curr = pt(p).ieeg;
    nfiles = length(curr.file);

    sleep_anns = {};
    sleep_times = [];
    sleep_files = [];

    wake_anns = {};
    wake_times = [];
    wake_files = [];

    % Loop over files
    for f = 1:nfiles

        % Loop over annotations
        nann = length(curr.file(f).ann);
        for a = 1:nann
            
            if strcmp(curr.file(f).ann,'empty'), continue; end
            
            % Loop over events
            nev = length(curr.file(f).ann(a).event);
            for e = 1:nev

                % Get the annotation
                annotation = curr.file(f).ann(a).event(e).description;
                time = curr.file(f).ann(a).event(e).start;

                % See if it contains text suggesting sleep (specifically N2 or
                % N3) vs wake

                %% Sleep
                if contains(annotation,'sleep','ignorecase',true) || ...
                        contains(annotation,'N2','ignorecase',true) || ...
                        contains(annotation,'N3','ignorecase',true) || ...
                        contains(annotation,'sws','ignorecase',true) || ...
                        contains(annotation,'spindle','ignorecase',true)
                    sleep_anns = [sleep_anns;annotation];
                    sleep_times = [sleep_times;time];
                    sleep_files = [sleep_files;f];

                %% Wake
                elseif contains(annotation,'awake','ignorecase',true) || ...
                        contains(annotation,'blink','ignorecase',true) || ...
                        contains(annotation,'pdr','ignorecase',true)
                    wake_anns = [wake_anns;annotation];
                    wake_times = [wake_times;time];
                    wake_files = [wake_files;f];

                end

            end
        end

    end

    %% Turn these into tables
    sleepT = table(sleep_files,sleep_times,sleep_anns);
    wakeT = table(wake_files,wake_times,wake_anns);

    %% Add to pt file
    pt(p).sw.sleep = sleepT;
    pt(p).sw.wake = wakeT;
end

end