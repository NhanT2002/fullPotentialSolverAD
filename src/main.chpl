use IO;
use Time;
use Math;
import input.potentialInputs;
use mesh;
use spatialDiscretization;
use temporalDiscretization;
use testMetrics;
use linearAlgebra;

config const runTests = false;  // Use --runTests=true to enable verification

proc main() {
    
    if runTests {
        run_tests();
    }
    else {
        writeln("Full Potential Solver - Chapel Implmentation");
        var t_ini: stopwatch;
        t_ini.start();

        var inputs = new potentialInputs();
        inputs.initializeFlowField();

        var Mesh = new shared MeshData(inputs.GRID_FILENAME_, inputs.ELEMENT_TYPE_);
        Mesh.buildConnectivity();
        writeln("Mesh loaded: ", Mesh.nelem_, " elements, ", Mesh.nedge_, " faces.");

        if inputs.FLOW_ == "steady" {
            var spatialDisc = new shared spatialDiscretization(Mesh, inputs);
            var steadySolver = new shared temporalDiscretization(spatialDisc.borrow(), inputs);
            
            steadySolver.initialize();
            steadySolver.solve();
        }
        else if inputs.FLOW_ == "IBM" {
            var spatialDisc = new shared spatialDiscretizationIBM(Mesh, inputs);

            var steadySolver = new shared temporalDiscretization(spatialDisc.borrow(), inputs);
            steadySolver.initialize();
            steadySolver.solve();
        }
        else {
            writeln(inputs.FLOW_, " flow type not recognized. Please check the configuration.");
        }

        t_ini.stop();
        writeln("total execution : ", t_ini.elapsed(), " seconds");
    }
    
}
