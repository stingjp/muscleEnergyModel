function [model] = probeActivate(model)
    import org.opensim.modeling.*
    
    % get the probeset
    probeset = model.getProbeSet();
    numProbes = probeset.getSize();

    % need to loop through and set them all to be enabled hopefully
    for p = 0:numProbes-1
        probe = probeset.get(p);
        probe.setEnabled(true);
    end

    % update the model to be returned. 
    model.updProbeSet();
end
