function runRRA(filename)
    import org.opensim.modeling.*
    
    rra = RRATool(filename);
    rra.run();

    disp('RRA run successful... files in results.')

end