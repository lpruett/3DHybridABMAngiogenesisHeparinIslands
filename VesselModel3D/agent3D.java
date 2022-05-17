package VesselModel3D;

import HAL.GridsAndAgents.SphericalAgent3D;
import HAL.Util;

/**
 * IN ORDER FOR PROGRAM TO WORK, MUST MODIFY SphericalAgent3D class.
 * See Required_additions.txt and the ReadMe for specifics.
 */

public class agent3D extends SphericalAgent3D<agent3D, grid3D> {
    ////////////////
    // PROPERTIES //
    ////////////////

    int color; // the color of the agent
    int current_persistency_time;  // time since the agent last chose new a direction of chemotaxis (ticks)
    double elongation_length; // the current length that the vessel has grown in it's current direction
    boolean heparinOn = false; // Whether heparin is releasing VEGF or not (true: releasing, false: not releasing)
    int timeSinceLastBranch; // The time since the vessel has last branched (ticks)
    boolean noOverlap = false; // determines if a body cell is overlapping with a MAP particle or not: if so it should move to not overlap (implementation not functioning, still in debugging process)
    double branching_probability; // probability that this vessel with branch (determined by VEGF presence)
    double[] lastPlace; // last place that a body cell was placed by a head cell
    String side = ""; // determines whether the cell originates from the "L" or "R" side (L side Z=0, R is side Z=Z)

    // Heparin MicroIslands (for making surface gradients)
    private int[] zero_VEGF; // Neighborhood inside the heparin island that should not have any gradient
    private int[] add_VEGF; // The neighborhood that includes the entirety of the particle

    // ^^^ Note that for the above neighborhoods, they will be used to make the Heparin MicroIslands only have gradient
    // on the edges.  The zero_VEGF neighborhood will be of radius slightly smaller than the add_VEGF neighborhood.
    // They will be concentric.  This means that they overlap for most of the center, but for a ring outside the zero_VEGF
    // neighborhood, there will be only add_VEGF.  The zero_VEGF will constantly remove VEGF from the center of the particle
    // after the add_VEGF continuously adds.  This means that there will only be a ring around the heparin MicroIsland
    // that actually has VEGF.


    ////////////////////
    // PROPERTY TYPES //
    ////////////////////

    // CELL TYPES
    public static final int HEAD_CELL = 1;
    public static final int BODY_CELL = 2;
    public static final int MAP_PARTICLE = 3;
    public static final int HEPARIN_ISLAND = 4;
    public static final int MACROPHAGE = 5;

    // COLORS
    public static int HEAD_CELL_COLOR = Util.RED;
    public static int BODY_CELL_COLOR = Util.RED;
    public static int MAP_PARTICLE_COLOR = Util.RGB(23.0/255, 28.0/255, 173.0/255); // normal MAP;
    public static int HEPARIN_MAP_COLOR = Util.RGB(48.0/255, 191.0/255, 217.0/255); // Heparin MAP;
    public static int MACROPHAGE_COLOR = Util.WHITE;

    ///////////////////////
    // UTILITY VARIABLES //
    ///////////////////////

    public static final double MAX_VELOCITY = grid3D.MIGRATION_RATE; // Max speed of the head cell as it migrates.
    public static final double SCALE_FACTOR = grid3D.SCALE_FACTOR;
    public static final double TIME_CALE_FACTOR = grid3D.TIME_SCALE_FACTOR;

    int[] posneg = {1, -1}; // utility variable that contains 1 and -1

    // SEE grid3D.java for variable descriptions
    public static double MAP_RAD = grid3D.MAP_RAD;
    public static final double VESSEL_RADIUS = grid3D.VESSEL_RADIUS;
    public static final double VESSEL_VEGF_CONSUME = grid3D.VESSEL_VEGF_CONSUME;
    public static final int  AGE_BEFORE_CONSUME = grid3D.AGE_BEFORE_CONSUME;
    public static final double MIGRATION_RATE = grid3D.MIGRATION_RATE;
    public static final double VEGF_SENSITIVITY_THRESHOLD = grid3D.VEGF_SENSITIVITY_THRESHOLD;
    public static final double MAX_ELONGATION_LENGTH = grid3D.MAX_ELONGATION_LENGTH;
    public static final double PERSISTENCY_TIME = grid3D.PERSISTENCY_TIME;
    public static final double BRANCH_DELAY = grid3D.BRANCH_DELAY;
    public final static double LOW_BRANCHING_PROBABILITY= grid3D.LOW_BRANCHING_PROBABILITY;
    public final static double LOW_MED_VEGF_THRESHOLD = grid3D.LOW_MED_VEGF_THRESHOLD;
    public final static double MED_BRANCHING_PROBABILITY= grid3D.MED_BRANCHING_PROBABILITY;
    public final static double MED_HIGH_VEGF_THRESHOLD = grid3D.MED_HIGH_VEGF_THRESHOLD;
    public final static double HIGH_BRANCHING_PROBABILITY= grid3D.HIGH_BRANCHING_PROBABILITY;
    public static final double SEQUENTIAL_TURN_ON = grid3D.SEQUENTIAL_TURN_ON;
    public static double VEGF_SCALE_FACTOR = grid3D.VEGF_SCALE_FACTOR;


    ////////////////////
    // INITIALIZATION //
    ////////////////////

    /**
     * General agent initializing method
     * @param type agent type
     * @param radius the radius of the particle
     */
    public void Init(int type, double radius){

        this.radius = radius;
        this.type = type;

        if (type == HEAD_CELL) {
            this.color = HEAD_CELL_COLOR;
        } else if (type == BODY_CELL){
            this.color = BODY_CELL_COLOR;
        } else if (type == MAP_PARTICLE){
            this.color = MAP_PARTICLE_COLOR;
        } else if (type == HEPARIN_ISLAND) {
            this.color = HEPARIN_MAP_COLOR;
        } else if (type == MACROPHAGE) {
            this.color = MACROPHAGE_COLOR;
        }
    }

    /**
     * Initializes head cell with a specific side of origin
     */
    public void Init_HEAD_CELL(String side){
        this.type = HEAD_CELL;
        this.color = HEAD_CELL_COLOR;
        this.radius = VESSEL_RADIUS;
        this.elongation_length = 0;
        this.current_persistency_time = -1;
        this.timeSinceLastBranch = 0;
        this.lastPlace = new double[]{Xpt(), Ypt(), Zpt()};
        this.side = side;
    }

    /**
     * Initializes body cell
     */
    public void Init_BODY_CELL(){
        this.type = BODY_CELL;
        this.color = BODY_CELL_COLOR;
        this.radius = VESSEL_RADIUS;
    }

    /**
     * Called by the grid class.  Recursively generates MAP particles with the designated ratio of Heparin microIslands.
     * Makes sure that they do not overlap and are spaced according to HCP.
     */
    public void Recursive_MAP_Generator(double heparin_island_percentage){
        assert G != null;
        // recursive function, makes 12 surrounding MAP particles if the coordinates are in bounds and unoccupied
        double x_distance = grid3D.MAP_GAP_CENTERS * Math.cos(Math.toRadians(60));
        double y_distance = grid3D.MAP_GAP_CENTERS * Math.sin(Math.toRadians(60));

        // Make 6 surrounding MAP_seed on same plane
        for (int i : posneg) {
            for (int j : posneg) {
                // calculate new coordinate for the next MAP gel particle
                double new_agent_x = Xpt() + (i * x_distance);
                double new_agent_y = Ypt() + (j * y_distance);

                // make sure that the new coordinate is inside bounds
                boolean x_in_range = (0 < new_agent_x) && (new_agent_x < grid3D.x);
                boolean y_in_range = (0 < new_agent_y) && (new_agent_y < grid3D.y);

                // make sure that there is nothing there
                boolean open_for_MAP = false;
                if (x_in_range && y_in_range) {
                    int occlusions = G.PopAt((int) new_agent_x, (int) new_agent_y, (int) Zpt());
                    if (occlusions == 0) {
                        open_for_MAP = true;
                    }
                }

                // If the coordinate is within the range of the grid, and it is unoccupied, then make a new agent there.
                if (x_in_range && y_in_range && open_for_MAP) {
                    agent3D new_agent;

                    // has possibility of becoming a heparin MicroIsland or just a normal MAP particle
                    if(G.rng.Double() < heparin_island_percentage){
                        new_agent = G.NewAgentPT(new_agent_x, new_agent_y, Zpt());
                        new_agent.Init(HEPARIN_ISLAND, MAP_RAD);
                        int[] zeroHood = Util.SphereHood(true, MAP_RAD);
                        int[] oneHood = Util.SphereHood(false, MAP_RAD+1);
                        new_agent.zero_VEGF = zeroHood;
                        new_agent.add_VEGF = oneHood; //Together zero and add VEGF hoods allow for surface release of VEGF (only required for heparin MicroIslands, not for normal MAP gel)
                    } else {
                        new_agent = G.NewAgentPT(new_agent_x, new_agent_y, Zpt());
                        new_agent.Init(MAP_PARTICLE, MAP_RAD);
                    }
                    new_agent.Recursive_MAP_Generator(heparin_island_percentage); // make MAP paricles around the one just made
                }
            }
        }

        // 6 surrounding on either side of xy plane.  The 3 "above" and 3 "below" in a triangular formation

        double Zsc_distance = grid3D.MAP_GAP_CENTERS * Math.cos(Math.toRadians(30)); // distance above/below xy plane
        double Ysc_distance = grid3D.MAP_GAP_CENTERS * Math.sin(Math.toRadians(30)); // distance above x coordinate

        // calculate position of MAP particle that is straight forward (.'.) quotation mark, if cells are r

        // calculate new coordinate
        double new_agent_y = Ypt() + (Ysc_distance);
        double new_agent_z = Zpt() + (Zsc_distance);

        // make sure that the new coordinate is inside bounds
        boolean y_in_range = (0 < new_agent_y) && (new_agent_y < grid3D.y);
        boolean z_in_range = (0 < new_agent_z) && (Zpt() < grid3D.z);


        // make sure that there is nothing there
        boolean open_for_MAP = false;
        if (y_in_range && z_in_range) {
            int occlusions = G.PopAt((int) Xpt(), (int) new_agent_y, (int) new_agent_z);
            if (occlusions == 0) {
                open_for_MAP = true;
            }
        }

        if (y_in_range && z_in_range && open_for_MAP){
            agent3D new_agent = G.NewAgentPT(Xpt(), new_agent_y, new_agent_z);
            new_agent.Init(MAP_PARTICLE, MAP_RAD);
            new_agent.Recursive_MAP_Generator(heparin_island_percentage);
        }
    }

    /////////////
    // METHODS //
    /////////////


    /**
     * Called on every agent by the grid class, specifying commands for each cell, depending on its type
     */
    public void StepCell() {
        assert G != null;

        // if a vessel reaches the center or a quarter, keep track of it for data export.
        if (this.side.equals("L")){
            if (Zpt() >= G.zDim/2.0 && grid3D.CenterArrivalTime == -1){
                grid3D.CenterArrivalTime = G.GetTick()/grid3D.TIME_SCALE_FACTOR;
            }
            if (Zpt() >= G.zDim/4.0 && grid3D.QuarterArrivalTime == -1){
                grid3D.QuarterArrivalTime = G.GetTick()/grid3D.TIME_SCALE_FACTOR;
            }
        } else if (this.side.equals("R")){
            if (Zpt() <= G.zDim/2.0 && grid3D.CenterArrivalTime == -1){
                grid3D.CenterArrivalTime = G.GetTick()/grid3D.TIME_SCALE_FACTOR;
            }
            if (Zpt() <= (G.zDim*3)/4.0 && grid3D.QuarterArrivalTime == -1){
                grid3D.QuarterArrivalTime = G.GetTick()/grid3D.TIME_SCALE_FACTOR;
            }
        }

        // HEAD CELLS
        if(type == HEAD_CELL) { // if type head cell

            // Consistently add to persistency time to keep track of how long this cell has persisted in this direction (-1 is a new head cell, so special case)
            if (current_persistency_time != -1){
                current_persistency_time += 1;
                timeSinceLastBranch += 1;
            }

            // Sensitivity Threshold: Stop everything if there is not enough VEGF (cell should not take any actions)
            if (G.VEGF.Get(Isq()) < VEGF_SENSITIVITY_THRESHOLD){
                return;
            }

            // if persistency time is -1, or both the persistency time and the elongation time has been reached... (persistency time = -1 is to notify newly initialized head cells to find a new direction)
            if ((current_persistency_time == -1)||((current_persistency_time > PERSISTENCY_TIME)&&(elongation_length+G.Dist(Xpt(), Ypt(), Zpt(), lastPlace[0], lastPlace[1], lastPlace[2]) >= MAX_ELONGATION_LENGTH))){
                // calculate branching probability
                CalculateBranchingProbability();

                // if it is time to branch,
                if (timeSinceLastBranch > BRANCH_DELAY){
                    if (G.rng.Double() < branching_probability) { // and if branch probability is satisfied, then
                        // make the vessel branch
                        double[] location = {Xpt() - 0.02, Ypt() - 0.02, Zpt() - 0.02};
                        if (G.In(location[0], location[1], location[2])) {
                            G.NewAgentPT(location[0], location[1], location[2]).Init_HEAD_CELL(this.side); // create another head cell
                            timeSinceLastBranch = 0; // change the branch time to 0
                        }
                    }
                }

                // find a new direction and update persistency times and elongation lengths.
                findNewDirection();
                current_persistency_time = 0;
                elongation_length = 0;
            }

            // make the head cell move forward towards it's target
            migrate_head();

        // BODY CELLS
        } else if (type == BODY_CELL){
            // if the body cell is overlapping with a MAP particle or a Heparin microIsland, then move!
            if (!noOverlap){
                double stillOverlap = SumForcesTyped(radius+MAP_RAD,(overlap, other)-> overlap, new int[] {MAP_PARTICLE, HEPARIN_ISLAND});
                ForceMove();
                if (stillOverlap == 0){ // If there is no longer overlap, then specify by changing noOverlap to true.  Will stop checking for body cell overlap with MAP particles.
                    noOverlap = true;
                }
            }

            if (this.Age() > AGE_BEFORE_CONSUME) {  // body cells that are of a certain tick age will consume VEGF (they wait before consuming so they don't disturb gradient around head cell)
                G.VEGF.Add(Isq(), -VESSEL_VEGF_CONSUME);
            }

        // HEPARIN MICROISLANDS
        } else if (type == HEPARIN_ISLAND){ // if the agent is a heparin microisland, then take appropriate actions (see stepHeparinIslands())
            stepHeparinIslands();
        }
    }

    /**
     * Finds the new gradient direction (called once head cell reaches elongation length and persistency time)
     */
    public void findNewDirection(){
        // GRADIENTS
        double gradX;
        double gradY;
        double gradZ;

        try{ // Make sure that the coordinate in question is still within bounds to prevent runtime errors
            assert G != null;
            gradX=G.VEGF.GradientX(Xpt(), Ypt(), Zpt());
            gradY=G.VEGF.GradientY(Xpt(), Ypt(), Zpt());
            gradZ=G.VEGF.GradientZ(Xpt(), Ypt(), Zpt());
        } catch (ArrayIndexOutOfBoundsException e) {
            return;
        }

        // CALCULATE MOVEMENT (velocity in direction of highest VEGF gradient)
        // Works by making a unit vector (divide each vector by the norm), and multiplying by the migration rate
        double norm= Util.Norm(gradX,gradY,gradZ);
        if(gradX!=0) {
            xVel = gradX / norm * MIGRATION_RATE;
        }
        if(gradY!=0) {
            yVel = gradY / norm * MIGRATION_RATE;
        }
        if(gradZ!=0) {
            zVel = gradZ / norm * MIGRATION_RATE;
        }
    }

    /**
     * Calls head cells to migrate in direction determined by findNewDirection
     */
    public void migrate_head(){

        double FORCE_SCALER = 1;

        double storeXVel = xVel;
        double storeYVel = yVel;
        double storeZVel  = zVel;

        // SUM FORCES AND MOVE (make sure that there is no overlap with MAP particle or Heparin microIslands)
        SumForcesTyped(radius+MAP_RAD,(overlap, other)-> overlap*FORCE_SCALER, new int[] {MAP_PARTICLE, HEPARIN_ISLAND});
        // TODO:  SumForcesTyped is not working as expected since correcting the migration rate from 3 to 30.  The vessels "bounce" off of the MAP particles and microislands instead of navigating around them.

        // MAKE VESSELS NOT CLUMP TOGETHER! stay away from each other.
        SumForcesTyped(10,(overlap, other)-> overlap*FORCE_SCALER, new int[] {HEAD_CELL});

        // max speed
        CapVelocity(MAX_VELOCITY);

        // call actual movement.
        ForceMove();

        xVel = storeXVel;
        yVel = storeYVel;
        zVel = storeZVel;

        // Leave a trail of body cells right behind (centers separated by distance=VESSEL_RADIUS)
        assert G != null;
        if (current_persistency_time != -1){
            if (G.Dist(Xpt(), Ypt(), Zpt(), lastPlace[0], lastPlace[1], lastPlace[2])>= VESSEL_RADIUS){
                if (G.In(Xpt()-.01, Ypt()-.01, Zpt()-.01)){
                    G.NewAgentPT(Xpt()-.01, Ypt()-.01, Zpt()-.01).Init_BODY_CELL();
                    lastPlace[0] = Xpt();
                    lastPlace[1] = Ypt();
                    lastPlace[2] = Zpt();
                    elongation_length += VESSEL_RADIUS;
                } else if (G.In(Xpt()+.01, Ypt()+.01, Zpt()+.01)){
                    G.NewAgentPT(Xpt()+.01, Ypt()+.01, Zpt()+.01).Init_BODY_CELL();
                    lastPlace[0] = Xpt();
                    lastPlace[1] = Ypt();
                    lastPlace[2] = Zpt();
                    elongation_length += VESSEL_RADIUS;
                }
            }
        }
    }

    /**
     * Calculates branching probability based on the VEGF present at the agent's coordinates.
     */
    public void CalculateBranchingProbability(){
        assert G != null;
        if (G.VEGF.Get(Isq()) < LOW_MED_VEGF_THRESHOLD){
            branching_probability = LOW_BRANCHING_PROBABILITY;
        } else if (G.VEGF.Get(Isq()) < MED_HIGH_VEGF_THRESHOLD){
            branching_probability = MED_BRANCHING_PROBABILITY;
        } else if (G.VEGF.Get(Isq()) > HIGH_BRANCHING_PROBABILITY){
            branching_probability = HIGH_BRANCHING_PROBABILITY;
        }
    }

    /**
     * Calls Heparin microIslands to release VEGF depending on presence of vessels nearby
     */
    public void stepHeparinIslands() {
        assert G != null;
        // (turn on means begin VEGF release)
        if (G.GetTick() < SEQUENTIAL_TURN_ON) { // if the first sequential turn on time has been reached, turn on the first eigth of microislands on each side
            if ((Zpt()<G.zDim/8.0)||(Zpt()>G.zDim*7/8.0)){
                heparinOn = true;
            }
        } else if (G.GetTick() < 2*SEQUENTIAL_TURN_ON){ // if the second sequential turn on time has been reached, turn on the first quarter of microislands on each side
            if ((Zpt()<G.zDim/4.0)||(Zpt()>G.zDim*3/4.0)){
                heparinOn = true;
            }
        } else if (G.GetTick() < 3*SEQUENTIAL_TURN_ON){ // if the third sequential turn on time has been reached, turn on the first 3/4 of microislands on each side
            if ((Zpt()<G.zDim*3/8.0)||(Zpt()>G.zDim*5/8.0)){
                heparinOn = true;
            }
        } else { // after the fourth sequential turn on, turn on all microislands
            heparinOn = true;
        }


        // release VEGF if necessary
        if (heparinOn){ // If heparin is supposed to release VEGF
            // sets edge of the Heparin islands to have gradient, but none on the inside. (Requires understanding of how HAL holds coordinates-- ask me for explanation!)

            // HAL holds coordinates in neighborhoods as [X1, Y1, Z1, X2, Y2, Z2, ... XN, YN, ZN]
            // neighborhood.length/3 gives how many coordinates in neighborhood
            // iterate through all points in neighborhood by iterating through length calculated above

            int len = add_VEGF.length/3;
            for (int i = 0; i < len; i++) {
                double[] location = {Xpt()+add_VEGF[3*i], Ypt()+add_VEGF[(3*i)+1], Zpt()+add_VEGF[(3*i)+2]};
                if (G.In(location[0], location[1], location[2])){
                    G.VEGF.Set(location[0], location[1], location[2], (1/grid3D.currentHeparinPercentage)*VEGF_SCALE_FACTOR);
                }
            }

            len = zero_VEGF.length/3;
            for (int i = 0; i < len; i++) {
                double[] location = {Xpt()+zero_VEGF[3*i], Ypt()+zero_VEGF[(3*i)+1], Zpt()+zero_VEGF[(3*i)+2]};
                if (G.In(location[0], location[1], location[2])){
                    G.VEGF.Set(location[0], location[1], location[2], 0);
                }
            }

            // If there is a head cell or body cell nearby, then stop releasing gradient
            for (agent3D agent : G.IterAgentsRad(Xpt(), Ypt(), Zpt(), MAP_RAD + VESSEL_RADIUS+1)) {
                if (agent.type == HEAD_CELL || agent.type == BODY_CELL) {
                    heparinOn = false;
                    break;
                }
            }
        }
    }
}


