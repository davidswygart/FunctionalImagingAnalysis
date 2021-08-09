function [attrs,epochs] = fastLoadSymphonyEpochAttributes(fname,protocolID)
%TODO: flesh this out more... in general we will want more metadata from
%symphony. See comments.

info = h5info(fname,'/');
epochs = [];

for eg = info.Groups.Groups(2).Groups' %epoch groups
    %eg.Attributes(...'label') == 'Control', etc.
    %eg.Links(...'source') OR eg.Groups(....'source')
    
    for eb = eg.Groups(1).Groups' %epoch blocks
        if contains(eb.Name,protocolID)
            epoch = eb.Groups(1).Groups(1);
            start_time = strcmp({epoch.Attributes(:).Name}, 'startTimeDotNetDateTimeOffsetTicks');
            end_time = strcmp({epoch.Attributes(:).Name}, 'endTimeDotNetDateTimeOffsetTicks');
            time_zone = strcmp({epoch.Attributes(:).Name}, 'startTimeDotNetDateTimeOffsetOffsetHours');

            st = datetime(uint64(epoch.Attributes(start_time).Value),'convertfrom','.net','timezone',num2str(epoch.Attributes(time_zone).Value));
            et = datetime(uint64(epoch.Attributes(end_time).Value),'convertfrom','.net','timezone',num2str(epoch.Attributes(time_zone).Value));

            tEpochs = repmat(cell2struct({epoch.Groups(2).Attributes(:).Value, st, et},{epoch.Groups(2).Attributes(:).Name, 't0', 'tf'},2),numel(eb.Groups(1).Groups),1);
            for e = 2:numel(eb.Groups(1).Groups)
                epoch = eb.Groups(1).Groups(e);
                start_time = strcmp({epoch.Attributes(:).Name}, 'startTimeDotNetDateTimeOffsetTicks');
                end_time = strcmp({epoch.Attributes(:).Name}, 'endTimeDotNetDateTimeOffsetTicks');
                time_zone = strcmp({epoch.Attributes(:).Name}, 'startTimeDotNetDateTimeOffsetOffsetHours');
                
                st = datetime(uint64(epoch.Attributes(start_time).Value),'convertfrom','.net','timezone',num2str(epoch.Attributes(time_zone).Value));
                et = datetime(uint64(epoch.Attributes(end_time).Value),'convertfrom','.net','timezone',num2str(epoch.Attributes(time_zone).Value));
                %epoch.Attributes : start time, end time
                %epoch.Groups(2).Attributes: epoch_params
                tEpochs(e) = cell2struct({epoch.Groups(2).Attributes(:).Value, st, et},{epoch.Groups(2).Attributes(:).Name, 't0', 'tf'},2);
            end
            for f = 1:numel(eb.Groups(2).Attributes)
                [tEpochs(:).(eb.Groups(2).Attributes(f).Name)] = deal(eb.Groups(2).Attributes(f).Value);
            end
            epochs = cat(1,epochs,tEpochs);
        end
    end
    
end
attrs = [];
% epochs = rmfield(epochs,{'symphonyVersion', 'protocolVersion','micronsPerPixel','angleOffsetFromRig'});
[~,i] = sort([epochs(:).t0],'ascend');
epochs = epochs(i);
ss = substruct('{}',{':'});
[epochs(:).t0] = subsref(num2cell(posixtime([epochs(:).t0])), ss);
[epochs(:).tf] = subsref(num2cell(posixtime([epochs(:).tf])), ss);

end