function spikeTimes = matchSpikes(epochs, cellDataPath)

c = load(cellDataPath).cellData;

dL = arrayfun(@(x) cell2mat(x.dataLinks.values()), c.epochs, 'uniformoutput',false);

spikeTimes = cellfun(@(x) spike_times(find(strcmp(dL,x))), {epochs(:).dataLink},'uniformoutput',false);


function y = spike_times(x)
if isempty(x)
    y = nan;
else
    y = c.epochs(x).get('spikes_ch1');
end
end

end
