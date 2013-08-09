package statalign.base;


import statalign.ui.McmcSettingsDlg;

/**
 * 
 * This is a container class containing the MCMC parameters described below.
 * 
 * @author miklos, novak
 *
 */
public class MCMCPars {
    public enum TreeType {
        STEINER("Steiner tree"),
        SPANNOID("Spannoid");

        private String label;
        private TreeType(String label) {
            this.label = label;
        }

        public String toString() {
            return label;
        }
    };

    /**
     * Tree topology type
     */
    public TreeType treeType;

    /**
     * Component size restriction when using Spannoids.
     */
    public int componentSize;

    /**
     * Restrict topology changes st. they conform to the component size restriction.
     */
    public boolean restrictTopologyChanges;

    /**
     * The strategy used in bonphy for edge contractions of Steiner Tree -> k-Spannoid.
     */
    public Spannoid.BonphyStrategy bonphyStrategy;

	/**
	 * 
	 * Number of burn-in steps
	 * 
	 */
	public int burnIn;
	/**
	 * Number of cycles total after burn-in period
	 * 
	 */
	public int cycles;       
	
	/**
	 * One sampling after each sampRate cycles
	 */
	public int sampRate;

	/**
	 * The seed for the random number generator
	 */
	public long seed;
	
	/**
	 * The swap seed for the swap random number generator.
	 */
	public long swapSeed;
	
	/**
	 * How often to propose swaps between chains (per cycle).
	 */
	public int swapRate;
	
	/**
	 * MCMC parameter automation settings
	 */
	public AutomateParamSettings autoParamSettings;
	
	/**
	 * This constructor sets the values in the class
	 * 
	 * @param burnIn this.burnIn is set to this value.
	 * @param cycles this.cycles is set to this value.
	 * @param sampRate this.sampRate is set to this value.
	 */
	public MCMCPars(TreeType treeType, int componentSize, boolean restrictTopologyChanges,
                    Spannoid.BonphyStrategy bonphyStrategy,
                    int burnIn, int cycles, int sampRate,
                    long seed, long swapSeed, int swapRate, AutomateParamSettings autoParamSettings) {
        this.treeType = treeType;
        this.componentSize = componentSize;
        this.restrictTopologyChanges = restrictTopologyChanges;
        this.bonphyStrategy = bonphyStrategy;
		this.burnIn = burnIn;
		this.cycles = cycles;
		this.sampRate = sampRate;
		this.swapRate = swapRate;
		this.seed = seed;
		this.swapSeed = swapSeed;
		this.autoParamSettings = autoParamSettings;
	}
	
	/*public MCMCPars(int burnIn, int cycles, int sampRate, long seed) {
		this.burnIn = burnIn;
		this.cycles = cycles;
		this.sampRate = sampRate;
		this.seed = seed;
	}*/
	
}
