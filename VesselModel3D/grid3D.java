package VesselModel3D;

import HAL.GridsAndAgents.AgentGrid3D;
import HAL.GridsAndAgents.PDEGrid3D;
import HAL.Gui.GridWindow;
import HAL.Gui.OpenGL3DWindow;
import HAL.Interfaces.DoubleToInt;
import HAL.Rand;
import HAL.Util;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.Timestamp;

public class grid3D extends AgentGrid3D<agent3D> {

    ////////////////
    // PARAMETERS //
    ////////////////

    // SCALE FACTORS
    public static final double SCALE_FACTOR = 0.1; // (distanceUnits/microns) microns to distance units
    public static final double TIME_SCALE_FACTOR = 10; // (ticks/hr) hours to ticks (normally 60?)
    public static final double VEGF_SCALE_FACTOR = 0.001; // unitless 0.0001

    // RUN SETTINGS
    public static final boolean EXPORT_DATA = true;
        public static final boolean EXPORT_HEAD_CELL_DATA = true; // requires EXPORT_DATA to be true
    public static final double HEAD_CELL_SAMPLE_TIME = 1 * TIME_SCALE_FACTOR; // (enter in hours) How frequently head cell distances from wound edges will be collected
    public static final double[] HEPARIN_ISLAND_PERCENTAGES = {0.1}; // enter as a decimal between 0 and 1, ratio of MAP particles that will be heparin microislands
    public static final double TRIALS = 1; // the number of trials that will be run per entry in HEPARIN_ISLAND_PERCENTAGES
    public final static int RUNTIME = (int)(168 *TIME_SCALE_FACTOR); // (enter time in hours) how long will the simulation run? - 168hr is 1 week




    // GRID PROPERTIES
    public static final int x = (int)(0.5 * (SCALE_FACTOR)*1000); // (enter distance in mm) dimension of the wound
    public static final int y = (int)(0.5 * (SCALE_FACTOR)*1000); // (enter distance in mm) dimension of the wound
    public static final int z = (int)(0.5 * (SCALE_FACTOR)*1000); // (enter distance in mm) dimension of the wound
    public static final double DIFFUSION_COEFFICIENT = 0.733; // (unitless) diffusion coefficient, ADI
    Rand rng = new Rand();

    // VIEW: what agents to display
    public static final boolean VIEW_MAP = false;  // Will standard MAP be visible?
    public static final boolean VIEW_HEP_ISLANDS = true; // Will heparin microislands be visible?
    public static final boolean VIEW_VESSELS = true; // Will vessels be visible?

    // MAP GEL
    public static final double MAP_DIAMETER = 80 * (SCALE_FACTOR); // (enter distance in microns)
    public static final double VESSEL_DIAMETER = 10 * (SCALE_FACTOR); // (enter distance in microns)
    public static final double MAP_GAP =  18 * (SCALE_FACTOR); // (enter distance in microns)
    public static final double SEQUENTIAL_TURN_ON = 24 * (TIME_SCALE_FACTOR); // (enter time in hours) time between each quadrant of heparin microislands being turned on

    // VESSELS
    public static final int NUM_VESSELS_PER_SIDE = 12; // The number of head vessels to start the model per side
    public static final double VESSEL_VEGF_CONSUME = 0.0001; // the amount of VEGF consumed by eligible cells (body cells older than AGE_BEFORE_CONSUME)
    public static final int AGE_BEFORE_CONSUME = 25; // (enter in ticks) age before a body cell can start consuming VEGF: to keep consumption from interacting with head cell gradient calculation
    public static final double MIGRATION_RATE = 30 * SCALE_FACTOR/TIME_SCALE_FACTOR; // (microns per hour) rate of vessel elongation
    public static final double VEGF_SENSITIVITY_THRESHOLD = 0.0001; // Threshold for VEGF sensitivity
    public static final double MAX_ELONGATION_LENGTH = 40 * (SCALE_FACTOR); // (enter in microns) max length a vessel can elongate in one direction before having to find another migration direction
    public static final double PERSISTENCY_TIME = 3 * (TIME_SCALE_FACTOR); // (enter in hours) the time a vessel must persist in the same direction before having to find another migration direction
    public static final double BRANCH_DELAY = 6 * (TIME_SCALE_FACTOR); // (enter in hours) The minimum amount hours between branching (model specific, included in Mehdizadeh et al.)
    // BRANCHING PROBABILITY AND THRESHOLDS_ PROBABILITIES NEED PARAMETERIZED BUT COULD STAY FIXED
    public final static double LOW_BRANCHING_PROBABILITY= 0.4; // probability of branching while VEGF is under LOW_MED_VEGF_THRESHOLD
    public final static double LOW_MED_VEGF_THRESHOLD = 0.05; // threshold of VEGF concentration between LOW and MED branching probabilities
    public final static double MED_BRANCHING_PROBABILITY= 0.55; // probability of branching while VEGF is between LOW_MED_VEGF_THRESHOLD and MED_HIGH_VEGF_THRESHOLD
    public final static double MED_HIGH_VEGF_THRESHOLD = 0.25; // threshold of VEGF concentration between MED and HIGH branching probabilities
    public final static double HIGH_BRANCHING_PROBABILITY= 0.9; // probability of branching while VEGF is above MED_HIGH_VEGF_THRESHOLD

    // DO NOT MODIFY
    public static final int HEAD_CELL = agent3D.HEAD_CELL;
    public static final int BODY_CELL = agent3D.BODY_CELL;
    public static final int MAP_PARTICLE = agent3D.MAP_PARTICLE;
    public static final int HEPARIN_ISLAND = agent3D.HEPARIN_ISLAND;

    public static final double MAP_RAD = (MAP_DIAMETER / 2.0);
    public static final double MAP_GAP_CENTERS = (MAP_GAP + (2 * MAP_RAD));
    public static final double VESSEL_RADIUS = VESSEL_DIAMETER/2.0;

    // DATA EXPORT
    public static StringBuilder CSV = new StringBuilder();  // contains the CSV that will contain general vessel metrics for data export
    public static StringBuilder HeadCellCSV = new StringBuilder(); // contains the distances of each head cell from their respective wound edge over time
    public static double CenterArrivalTime = -1;  // notes the time it takes for a single vessel to reach the center of the wound
    public static double QuarterArrivalTime = -1; // notes the time it takes for a single vessel to reach the first quarter of the wound

    // UTILITY
    public static double currentHeparinPercentage = 0;  // current heparin percentage (used in data export for file naming)
    public static int BatchNum = -1; // current batch number (used in data export for file naming)
    public static String date = "";  // current date (used in data export for file naming)
    public static StringBuilder percentages = new StringBuilder();  // list of heparin percentages in string form (used in data export for file naming)
    /////////////////
    // MAIN METHOD //
    /////////////////

    public static void main(String[] args) throws IOException{

        String folderName = "";
        if (EXPORT_DATA){  // if Data export is set to true, then prepare CSV and folders for data export
            // Verify that folder exists, else make a folder
            folderName = MakeFolder();
            // INITIALIZE CSV
            Initialize_CSV();
        }

        // INITIALIZE WINDOWS (used for visualization of vessels and VEGF)
        OpenGL3DWindow window = new OpenGL3DWindow("Angiogenesis", 900, 900, x, y, z);
        GridWindow VEGF_xz = new GridWindow("VEGF Diffusion X-Z plane", x, z,3);

        // INITIALIZE GRID THAT HOUSES AGENTS
        grid3D woundGrid = new grid3D(x, y, z);

        // RUN SIMULATIONS
        for (double heparinIslandPercentage : HEPARIN_ISLAND_PERCENTAGES) {  // for each percentage
            currentHeparinPercentage = heparinIslandPercentage; // remember what percentage is currently looked at for data export
            for (int trial = 0; trial < TRIALS; trial++) { // for each trial

                // RESET MODEL before starting a new trial
                woundGrid.VEGF.SetAll(0); // make sure that the VEGF grid is cleared
                woundGrid.ResetHard(); // reset the grid3D to remove all agents

                int one_indexed_trial = trial+1;  // give a 1 indexed value for trials (e.g. 0 becomes 1, 1 becomes 2, etc.)

                System.out.println("Started Trial " + one_indexed_trial + " of " + heparinIslandPercentage*100 + "% heparin islands..."); // system response that confirms model has started running

                System.out.println("...");

                // if data export is True, then initialize Head Cell data CSV
                if (EXPORT_DATA && EXPORT_HEAD_CELL_DATA){
                    Initialize_Head_CSV(one_indexed_trial, heparinIslandPercentage);
                }


                // INITIALIZE MAP PARTICLES AND VESSELS
                Init_MAP_Particles(woundGrid, heparinIslandPercentage); // initializes vessels in HCP (hexagonal compact packing)
                Init_Vessels(woundGrid);

                // TICK ACTIONS
                for (int step = 0; step < RUNTIME; step++) {  // For each timestep
                    woundGrid.StepVEGF(); // Step the gradient (allowing it to diffuse)
                    woundGrid.StepCells(); // Step all the cells
                    woundGrid.DrawGrid(window); // Draw the updated window
                    woundGrid.DrawGradientWindowed(VEGF_xz, Util::HeatMapBGR); // Draw the updated PDE grid
                    woundGrid.DrawAgents(window); // draw the new updated agents
                    woundGrid.VEGF.Update(); // Update the VEGF window with the newly drawn PDE grid
                    window.Update();// Update the wound grid window with the newly drawn agents
                    woundGrid.HeadCellDistOverTime(); // sample the head distance data if necessary

                    woundGrid.IncTick(); // Increment the time tick

                    if(window.IsClosed()){ // If the window is X'ed out, then
                        window.Close(); // close the wound grid window,
                        VEGF_xz.Close(); // and close the PDE window
                        break; // exit the time tick loop and end the program
                    }
                }

                System.out.println("Completed Trial " + one_indexed_trial + " of " + heparinIslandPercentage*100 + "% heparin islands."); // system response once each trial is completed

                // COLLECT DATA
                CollectVesselData(woundGrid, heparinIslandPercentage, one_indexed_trial);
            }
        }

        // if EXPORT_DATA is true, then export data (same with EXPORT_HEAD_CELL_DATA)
        if (EXPORT_DATA){
            ExportData(folderName);
            if (EXPORT_HEAD_CELL_DATA){
                ExportHeadCellTimeData(folderName);
            }
        }

        // Close the window
        window.Close(); // close the wound grid window
        VEGF_xz.Close(); // and close the PDE window
    }

    ///////////////
    // VARIABLES //
    ///////////////

    PDEGrid3D VEGF; // The PDE grid used to model VEGF diffusion

    ///////////////////
    // NEIGHBORHOODS //
    ///////////////////



    //////////////////
    // CONSTRUCTORS //
    //////////////////

    /**
     * Constructs the grid
     * @param x The x dimension of the wound
     * @param y The y dimension of the wound
     * @param z The z dimension of the wound
     */
    public grid3D (int x, int y, int z) {
        super(x, y, z, agent3D.class);
        VEGF = new PDEGrid3D(x, y, z);

    }

    ////////////////////
    // INITIALIZATION //
    ////////////////////

    /**
     * Initializes the wound grid with MAP particles by starting with a randomly placed seed MAP, and recursively generating MAP
     * around it in HCP packing
     * @param grid
     */
    public static void Init_MAP_Particles(grid3D grid, double heparin_island_percentage){
        agent3D MAP_seed = grid.NewAgentPTSafe(x * Math.random(), y * Math.random(), 0);  // Generates the "seed" particle
        MAP_seed.Init(MAP_PARTICLE, MAP_RAD); // Initializes the seed as a MAP particle
        MAP_seed.Recursive_MAP_Generator(heparin_island_percentage); // Creates all the particles around the seed, recursively (defined in agent3D class)
    }

    /**
     * Initializes the host vasculature in a given grid by placing head cells on opposite wound edges
     * @param grid The grid that the vasculature is to be initialized in
     */
    public static void Init_Vessels(grid3D grid){

        // first side (L)
        boolean empty = true;  // is the location occupied with MAP gel or not? (assume is empty=true)
        for (int i = 0; i < NUM_VESSELS_PER_SIDE;) {  // for as many vessels you want to start with  ("i" is tally for how many successful head vessels have been initialized),
            empty = true; // assume that the desired location is empty
            double[] location = {(x/2.0)*grid.rng.Double()+(x/4.0), (y/2.0)*grid.rng.Double()+(y/4.0), 1.5*VESSEL_RADIUS}; // starts at a random point on the z=1.5*vessel_radius plane (i.e the beginning of the wound, but with some leeway)
            for (agent3D agent : grid.IterAgentsRad(location[0], location[1], location[2], MAP_RAD+VESSEL_RADIUS)) { // Iterate through all locations around the desired point in a radius equal to MAP radius
                if (agent.type == MAP_PARTICLE || agent.type == HEPARIN_ISLAND) { // If there is a MAP particle center in that radius (meaning that it overlaps with the desired location), then
                    empty = false; // that means that the desired location is not empty. So,
                    break; // exit the for loop since you know that it is not empty
                }
            }
            if (empty){ // BUT if it is empty, initialize a head vessel there and increment "i" which is a tally for how many vessels have been initialized
                grid.NewAgentPT(location[0], location[1], location[2]).Init_HEAD_CELL("L"); // initialize a head cell at this location, and
                i++; // increment i.
            }
        }

        // other side (R)
        for (int i = 0; i < NUM_VESSELS_PER_SIDE;) {  // for as many vessels you want to start with  ("i" is tally for how many successful head vessels have been initialized)
            empty = true; // assume that the desired location is empty
            double[] location = {(x/2.0)*grid.rng.Double()+(x/4.0), (y/2.0)*grid.rng.Double()+(y/4.0), z-(1.5*VESSEL_RADIUS)}; // starts at a random point on the z=1.5*vessel_radius plane (i.e the beginning of the wound, but with some leeway)
            for (agent3D agent : grid.IterAgentsRad(location[0], location[1], location[2], MAP_RAD+VESSEL_RADIUS)) { // Iterate through all locations around the desired point in a radius equal to MAP radius
                if (agent.type == MAP_PARTICLE || agent.type == HEPARIN_ISLAND) { // If there is a MAP particle center in that radius (meaning that it overlaps with the desired location), then
                    empty = false; // the desired location is not empty. So,
                    break; // exit the for loop since you know that it is not empty
                }
            }
            if (empty){ // BUT if it is empty, initialize a head vessel there and increment "i" which is a tally for how many vessels have been initialized
                grid.NewAgentPT(location[0], location[1], location[2]).Init_HEAD_CELL("R"); // initialize a head cell at this location, and
                i++; // increment i.
            }
        }
    }

    //////////////////
    // GRID ACTIONS //
    //////////////////

    /**
     * Called on all cells inside the grid: calls StepCell on all cells, making each cell take their actions
     */
    public void StepCells() {
        for (agent3D cell: this) { // for each of the cells in the grid
            cell.StepCell(); // call them to take action
        }
    }

    /**
     * Calls the VEGF grid to take action, i.e. diffusing by one "tick"
     */
    public void StepVEGF(){
        VEGF.DiffusionADI(DIFFUSION_COEFFICIENT);
    }

    //////////////////
    // GRID DRAWING //
    //////////////////

    /**
     * Draws the grid inside the given window (just the background and outline of the used space)
     * @param window The window in which to draw the grid
     */
    public void DrawGrid (OpenGL3DWindow window){
        window.ClearBox(Util.RGB(225/245.0,198/245.0,153/245.0), Util.BLUE); // draw the grid with background color and line color
    }

    /**
     * Draw Agents depending on their color as spheres
     * @param window The window to draw the agents in
     */
    public void DrawAgents(OpenGL3DWindow window){
        for (agent3D cell : this) {
            if ((cell.type == MAP_PARTICLE) && (!VIEW_MAP)){
                continue;
            } else if ((cell.type == HEPARIN_ISLAND) && (!VIEW_HEP_ISLANDS)) {
                continue;
            } else if (((cell.type == HEAD_CELL) || (cell.type == BODY_CELL)) && (!VIEW_VESSELS)){
                continue;
            }
            window.CelSphere(cell.Xpt(),cell.Ypt(),cell.Zpt(),cell.radius, cell.color);
        }
    }

    /**
     * Draws the gradient in a separate window from the wound grid (gradient visualized from a top-down view (x-z plane))
     * by compressing the y dimension.
     * @param window The window to draw the gradient in
     * @param DrawConcs Use Util.HeatMap___ as the argument to visualize the diffusion
     */
    public void DrawGradientWindowed (GridWindow window, DoubleToInt DrawConcs){
        for (int x = 0; x < VEGF.xDim; x++) {
            for (int z = 0; z < VEGF.zDim; z++) {
                double VEGF_Sum=0;
                //add column to avgConcs
                for (int y = 0; y < VEGF.yDim; y++) { // add all the concentrations in a column (all y coordinates for intersection of the x-z plane)
                    VEGF_Sum+=VEGF.Get(x,y,z);
                }
                window.SetPix(x,z,DrawConcs.DoubleToInt(VEGF_Sum)); // draw the concentration to the x-y grid
            }
        }
    }

    /**
     * Draws the diffusion gradient on the bottom face of the wound grid visualization (gradient visualized from a top-down view (x-z plane))
     * by compressing the y dimension.
     * @param window The window containing the wound grid
     * @param DrawConcs Use Util.HeatMap___ as the argument to visualize the diffusion
     */
    public void DrawGradientEmbedded (OpenGL3DWindow window, DoubleToInt DrawConcs){
        window.ClearBox(Util.RGB(225/245.0,198/245.0,153/245.0), Util.BLUE);
        for (int x = 0; x < VEGF.xDim; x++) {
            for (int z = 0; z < VEGF.zDim; z++) {
                double VEGF_Sum=0;
                //add column to avgConcs
                for (int y = 0; y < VEGF.yDim; y++) { // add all the concentrations in a column (all y coordinates for intersection of the x-z plane)
                    VEGF_Sum+=VEGF.Get(x,y,z);
                }
                VEGF_Sum/=VEGF.yDim;
                window.SetPixXZ(x,z,DrawConcs.DoubleToInt(VEGF_Sum)); // draw the concentration to the bottom face of the model
            }
        }
    }

    /////////////////
    // DATA EXPORT //
    /////////////////

    /**
     * Initializes the main CSV (contained in variable, CSV) with data headers
     */
    public static void Initialize_CSV(){
        CSV.append("Trial Number, Heparin Percentage (%), Total BVL (microns), Inner Quadrant BV percentage, Outer Quadrant BV percentage, Quarter Arrival Time (h), Center Arrival Time (h)");
    }


    /**
     * Adds data headers the the Head Cell Data CSV file
     * @param trial the current trial number
     * @param heparinIslandPercentage the current heparin island percentage
     */
    public static void Initialize_Head_CSV(int trial, double heparinIslandPercentage) {
        if (HeadCellCSV.length() == 0){ // if the CSV is empty, then just list the data header information
            HeadCellCSV.append("TRIAL ").append(trial).append(", Heparin Island Percentage: ").append(Math.round(heparinIslandPercentage * 100)).append("\nTime (h), Head Cell Distances from wound edge (microns)");
        } else { // if the CSV is not empty, then make new lines to separate the trials, then list the data header information
            HeadCellCSV.append("\n\n\nTRIAL ").append(trial).append(", Heparin Island Percentage: ").append(Math.round(heparinIslandPercentage * 100)).append("\nTime (h), Head Cell Distances from wound edge (microns)");
        }
    }


    /**
     * Checks if the Data Export folder exists.  If it doesn't, then it makes one.
     * @return the file location of the folder that data should be exported to
     * @throws IOException if a folder is not able to be made.
     */
    public static String MakeFolder() throws IOException {
        Path folderName= Path.of("VesselModel3D\\Model3D_Data");  // holds the path that should lead to the data export folder
        File Model3D_Datafile = new File(String.valueOf(folderName));
        if (!Model3D_Datafile.exists()) { // if this folder does not exist, then
            if (!Model3D_Datafile.mkdir()) { // make the folder.  And if the folder cannot be made, then
                throw new IOException("Model3D_Data folder not made"); // throw an error.
            }
        }

        // Make a folder for all data taken on this date
        Timestamp timestamp = new Timestamp(System.currentTimeMillis());  // keeps track of the current timestamp
        String readableTime = ((timestamp.toString().replace(" ","_").replace(".", "-").replace(":", "-")).substring(0, 10)); // convert the timestamp to a readable year-month-date
        date = readableTime; // hold on to the current time stamp for future use.
        String timestamp_string = "VesselModel3D\\Model3D_Data\\" + readableTime; // holds the path that should lead to all data collected on this given day
        Path dateFolderName= Path.of(timestamp_string);
        File dateDataFolder = new File(String.valueOf(dateFolderName));
        if (!dateDataFolder.exists()) { // if this folder does not exist, then
            if (!dateDataFolder.mkdir()) {  // make the folder.  And fi the folder cannot be made, then
                throw new IOException("Date folder not made"); // throw an error.
            }
        }

        // convert the Heparin percentages into a readable string, then keep the value in "percentages" for future use.
        percentages.append("[");
        for (double percentage : HEPARIN_ISLAND_PERCENTAGES) {
            percentage = Math.round(percentage*100);
            if (percentages.length() != 1) {
                percentages.append(", ");
            }
            percentages.append((int)percentage);
        }
        percentages.append("]");


        // Make a batch folder that contains all data for this batch. (follows same logic as past, except labels batches with 1, 2, 3, etc. depending on what already exists.)
        int batchNum = 1;
        String batchString = timestamp_string + "\\"+ percentages +" Batch " + batchNum;
        Path batchFolderName= Path.of(batchString);
        while (Files.exists(batchFolderName)){
            batchNum ++;
            batchString = timestamp_string + "\\"+ percentages +" Batch " + batchNum;
            batchFolderName= Path.of(batchString);
        }
        BatchNum = batchNum;
        File batchDataFolder = new File(String.valueOf(batchFolderName));
        if (!batchDataFolder.mkdir()) {
            throw new IOException("Batch folder not made");
        }
        return batchString;
    }


    /**
     * Collect data on percent vessels in inner 50% and outer 25% on each side, quarter arrival time, and center arrival time.
     * @param G the grid that the vessels and MAP gel exist in
     * @param heparinIslandPercentage the heparin island percentage of the current trial
     * @param trial the current trial number
     */
    public static void CollectVesselData(grid3D G, double heparinIslandPercentage, int trial){

        // check their quadrants
        int outerQuadrant = 0; // counter for agents in the outer half of the wound (0 to x/2 since x is wound center)
        int innerQuadrant = 0; // counter for agents in the inner half of the wound (x/2 to x)
        int numHeadCells = 0; // counter for number of head cells


        // iterate through all agents, and keep track of the number of head cells, and total number of head/body cells.
        // also keep track how many cells re in the inner/outer quadrants.
        for (agent3D agent3D : G.IterAgentsRect(0, 0, 0, x, y, z)) {
            if (agent3D.type == HEAD_CELL){
                numHeadCells ++;
            }
            if ((agent3D.type == HEAD_CELL) || (agent3D.type == BODY_CELL)){
                if (agent3D.Zpt() < 3*z/4.0 && agent3D.Zpt() > z/4.0) {
                    innerQuadrant ++;
                } else{
                    outerQuadrant++;
                }
            }
        }
        int totalAgents = outerQuadrant + innerQuadrant;

        double ratioInner = innerQuadrant/(totalAgents*1.0);
        double ratioOuter = outerQuadrant/(totalAgents*1.0);
        double totalVesselLength = ((totalAgents*VESSEL_RADIUS) + (numHeadCells*VESSEL_RADIUS))/SCALE_FACTOR;  // since cells overlap by 50%, each body cell only contributes VESSEL_RADIUS to the vessel length
                                                                                                                // also, since head cells contribute the whole vessel diameter, add another 50% for each head cell

        // append the data to the main CSV file
        CSV.append("\nTrial ").append(trial).append(",").append(heparinIslandPercentage).append(", ").append(totalVesselLength).append(", ").append(ratioInner).append(", ").append(ratioOuter).append(", ").append(QuarterArrivalTime).append(", ").append(CenterArrivalTime);

        // Reset Arrival Times
        QuarterArrivalTime = -1;
        CenterArrivalTime = -1;

    }


    /**
     * Collect data on head cell distances from wound edges at appropriate time ticks (as dictated by HEAD_CELL_SAMPLE_TIME)
     */
    public void HeadCellDistOverTime(){
        // if head cell data is being exported
        if (EXPORT_HEAD_CELL_DATA){
            if (GetTick()% HEAD_CELL_SAMPLE_TIME == 0){ // and the sample time has been reached
                HeadCellCSV.append("\n").append(GetTick()/TIME_SCALE_FACTOR); // get the current time in hours and append it to the CSV

                // look through all cells
                for (agent3D agent3D : this.IterAgentsRect(0, 0, 0, x, y, z)) {
                    // get each head cell's distance from their respective sides, and append them to the CSV.
                    if (agent3D.type == HEAD_CELL){
                        if (agent3D.side.equals("L")){
                            HeadCellCSV.append(",").append(agent3D.Zpt()/(SCALE_FACTOR));
                        } else if (agent3D.side.equals("R")) {
                            HeadCellCSV.append(",").append((z-agent3D.Zpt())/(SCALE_FACTOR));
                        }
                    }
                }
            }
        }
    }


    /**
     * Export data from the main CSV into the given folder
     * @param folderName the name of the folder that the datafile is to be placed in
     * @throws IOException if the file cannot be written
     */
    public static void ExportData(String folderName) throws IOException {
        Path fileName3D= Path.of(folderName + "\\OverallBatchData(B"+BatchNum+"_"+date+"_" + percentages +").csv");
        Files.writeString(fileName3D, CSV);
    }

    /**
     * Export data from the head cell data CSV into the given folder
     * @param folderName the name of the folder that the datafile is to be placed in
     * @throws IOException if the file cannot be written
     */
    public static void ExportHeadCellTimeData (String folderName) throws IOException {
        Path fileNameHead= Path.of(folderName + "\\HeadCellTimeData(B"+BatchNum+"_"+date+"_" + percentages +").csv");
        Files.writeString(fileNameHead, HeadCellCSV);
    }
}