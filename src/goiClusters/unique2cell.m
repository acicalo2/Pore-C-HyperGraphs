function unique_genes = unique2cell(genes)
    [r,c]=size(genes);
    genes=reshape(genes,[r*c,1]);
    idx=cellfun(@(x) ~isempty(x),genes);
    unique_genes=unique(genes(idx));
end