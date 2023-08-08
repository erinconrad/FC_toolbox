function out = convert_rid_hupid(id,which,conversion_table)

switch which
    case 'rid'
        row = conversion_table.record_id == id;
        out = conversion_table.hupsubjno(row);
    case 'hupid'
        row = conversion_table.hupsubjno == id;
        out = conversion_table.record_id(row);

end

end