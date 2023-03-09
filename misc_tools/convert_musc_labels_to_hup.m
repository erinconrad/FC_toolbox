function hup_labels = convert_musc_labels_to_hup(musc_labels)

hup_labels = strrep(musc_labels,'AH','B');
hup_labels = strrep(hup_labels,'PH','C');


end