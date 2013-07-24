package statalign.base;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Locale;
import java.util.Random;

import mpi.MPI;
import statalign.MPIUtils;
import statalign.base.thread.Stoppable;
import statalign.base.thread.StoppedException;
import statalign.distance.Distance;
import statalign.postprocess.PostprocessManager;
import statalign.postprocess.plugins.contree.CNetwork;
import statalign.ui.ErrorMessage;
import statalign.ui.MainFrame;
import statalign.utils.SimpleStats;

import com.ppfold.algo.AlignmentData;
import com.ppfold.algo.FuzzyAlignment;

/**
 * 
 * This class handles an MCMC run.
 * 
 * The class extends <tt>Stoppable</tt>, it may be terminated/suspended in
 * graphical mode.
 * 
 * @author miklos, novak
 * 
 */
public class Mcmc extends Stoppable {

	// Constants

	// int samplingMethod = 1; //0: random sampling, 1: total sampling
	final static int FIVECHOOSE[] = { 35, 5, 15, 35, 10 }; // edge, topology,
	// indel parameter,
	// alignment,
	// substitutionparameter
	final static int FOURCHOOSE[] = { 35, 5, 25, 35 }; // edge, topology, indel
	// parameter, alignment
	private static final int SAMPLE_RATE_WHEN_DETERMINING_THE_SPACE = 100;
	private static final int BURNIN_TO_CALCULATE_THE_SPACE = 25000;


	// Parallelization

	/** Is this a parallel chain? By-default false. */
	private boolean isParallel = false;

	/** The number of processes. */
	private int noOfProcesses;

	/** The rank of the process. */
	private int rank;

	/** When this variable reaches 0 we do a swap. */
	private int swapCounter;

	/** The random number generator used for swapping. */
	private Random swapGenerator;

	// Non parallelization

	public CNetwork network; 

	/** Current tree in the MCMC chain. */
	// public ITree tree;

    private MCMCStrategy strategy;

	/**
	 * MCMC parameters including the number of burn-in steps, the total number
	 * of steps in the MCMC and the sampling rate.
	 */
	public MCMCPars mcmcpars;
	
	public AutomateParamSettings autoPar;

	public McmcStep mcmcStep = new McmcStep();

	/** PostprocessManager that handles the postprocessing modules. */
	public PostprocessManager postprocMan;

	/** True while the MCMC is in the burn-in phase. */
	public boolean burnin;

    private double heat;

	public Mcmc(MCMCStrategy strategy, MCMCPars mcmcpars, PostprocessManager ppm) {
		postprocMan = ppm;
		ppm.mcmc = this;
		this.strategy = strategy;
		this.mcmcpars = mcmcpars;
		this.autoPar = mcmcpars.autoParamSettings;
		heat = 1.0d;
	}

	public Mcmc(MCMCStrategy strategy, MCMCPars mcmcpars, PostprocessManager ppm,
			int noOfProcesses, int rank, double heat) {
		this(strategy, mcmcpars, ppm);
		this.noOfProcesses = noOfProcesses;
		this.rank = rank;
		heat = heat;

		// Is parallel!
		isParallel = true;
	}

	public int alignmentSampled = 0;
    public int alignmentAccepted = 0;
    public int edgeSampled = 0;
    public int edgeAccepted = 0;
    public int topologySampled = 0;
    public int topologyAccepted = 0;
    public int indelSampled = 0;
    public int indelAccepted = 0;
    public int substSampled = 0;
    public int substAccepted = 0;
	
	private static final DecimalFormat df = new DecimalFormat("0.0000");

	/**
	 * In effect starts an MCMC run. It first performs a prescribed number of
	 * burn-in steps, then, if one wants to automate the sampling rate, goes to 
	 * secondary burn-in where the sampling rate is determined. After that it
	 *  makes the prescribed number of steps after both burn-ins,
	 * drawing samples with the prescribes frequency. It also calls the
	 * appropriate functions of the PostpocessManager <tt>postprocMan</tt> to
	 * trigger data transfer to postprocessing modules when necessary
	 */
	public int doMCMC() {
		if (isParallel) {
			String str = String.format(
					"Starting MCMC chain no. %d/%d (heat: %.2f)\n\n", 
					rank + 1, noOfProcesses, heat);
			MPIUtils.println(rank, str);
			swapGenerator = new Random(mcmcpars.swapSeed);
		} else {
			System.out.println("Starting MCMC...\n");
		}

		MainFrame frame = postprocMan.mainManager.frame;
		Utils.generator = new Random(mcmcpars.seed + rank);
		long currentTime, start = System.currentTimeMillis();

		// Triggers a /before first sample/ of the plugins.
		if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
			postprocMan.beforeFirstSample();
		}

		ArrayList<Double> logLikeList = new ArrayList<Double>();
		int errorCode = 0;

		try {
			//only to use if AutomateParameters.shouldAutomate() == true
			ArrayList<String[]> alignmentsFromSamples = new ArrayList<String[]>(); 
			int burnIn = mcmcpars.burnIn;
			boolean stopBurnIn = false;


			if(autoPar.automateBurnIn){
				burnIn = 10000000;
			} 

			if(autoPar.automateSamplingRate){
				burnIn += BURNIN_TO_CALCULATE_THE_SPACE;
			}

			burnin = true;
			for (int i = 0; i < burnIn; i++) {

				sample(0);

				// Triggers a /new step/ and a /new peek/ (if appropriate) of
				// the plugins.
				if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
					// TODO do above inside sample() and add more info
					mcmcStep.newLogLike = getTree().getLogLike();
					mcmcStep.burnIn = burnin;
					postprocMan.newStep(mcmcStep);
					if (i % mcmcpars.sampRate == 0) {
						postprocMan.newPeek();
					}
				}
				
				//every 50 steps, add the current loglikelihood to a list
				// and check if we find a major decline in that list 
				if(autoPar.automateBurnIn && i % 50 == 0){
					logLikeList.add(getState().logLike);
					if(!stopBurnIn){
						stopBurnIn = AutomateParameters.shouldStopBurnIn(logLikeList);
						if(autoPar.automateSamplingRate && stopBurnIn){
							burnIn = i + BURNIN_TO_CALCULATE_THE_SPACE;
						}else if (stopBurnIn){
							burnIn = i;
						}
					}
				}
				currentTime = System.currentTimeMillis();
				int realBurnIn = burnIn - BURNIN_TO_CALCULATE_THE_SPACE;
				if (frame != null) {
					String text = "";
					if((i > realBurnIn ) && autoPar.automateSamplingRate){
						text = " Burn-in to aid automation of MCMC parameters: " + (i-realBurnIn + 1) ;
					}else{
						text = " Burn-in: " + (i + 1);
						if(autoPar.automateBurnIn)
							text += " -- total length will be determined automatically";
						else
							text += " out of "+burnIn;
					}
					frame.statusText.setText(text);
				} else if (i % 1000 == 999) {
					System.out.println("Burn in: " + (i + 1));
				}


				if( autoPar.automateSamplingRate && (i >= realBurnIn) && i % SAMPLE_RATE_WHEN_DETERMINING_THE_SPACE == 0)   {
					String[] align = getState().getLeafAlign();
					alignmentsFromSamples.add(align);
				}	
			}
			
			//both real burn-in and the one to determine the sampling rate have now been completed.
			burnin = false;

			int period;
			if(autoPar.automateNumberOfSamplesToTake){
				period = 1000000;
			}else{
				period = mcmcpars.cycles / mcmcpars.sampRate;
			}

			int sampRate;
			if(autoPar.automateSamplingRate){
				if(frame != null)
				{
					frame.statusText.setText(" Calculating the sample rate");
				}
				else
				{
					System.out.println(" Calculating the sample rate");
				}
				ArrayList<Double> theSpace = Distance.spaceAMA(alignmentsFromSamples);
				sampRate = AutomateParameters.getSampleRateOfTheSpace(theSpace,SAMPLE_RATE_WHEN_DETERMINING_THE_SPACE);

			}else{
				sampRate = mcmcpars.sampRate;
			}


			int swapNo = 0; // TODO: delete?
			swapCounter = mcmcpars.swapRate;
			AlignmentData alignment = new AlignmentData(getState().getLeafAlign());
			ArrayList<AlignmentData> allAlignments = new ArrayList<AlignmentData>();
			ArrayList<Double> distances = new ArrayList<Double>();

			boolean shouldStop = false;
			double currScore = 0;
			for (int i = 0; i < period && !shouldStop; i++) {
				for (int j = 0; j < sampRate; j++) {
					// Samples.
					sample(0);

					//FuzzyAlignment fuzzyAlignment2 = FuzzyAlignment.getFuzzyAlignmentAndProject(alignments, "");

					// Proposes a swap.
					if (isParallel) {
						swapCounter--;
						if (swapCounter == 0) {
							swapNo++;
							swapCounter = mcmcpars.swapRate;

							doSwap(swapNo);
						}
					}

					// Triggers a /new step/ and a /new peek/ (if appropriate)
					// of the plugins.
					if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
						mcmcStep.newLogLike = getTree().getLogLike();
						mcmcStep.burnIn = burnin;
						postprocMan.newStep(mcmcStep);
						if (burnIn + i * period + j % mcmcpars.sampRate == 0) {
							postprocMan.newPeek();
						}
					}

				}
				currentTime = System.currentTimeMillis();
				if (frame == null && !isParallel) {
					System.out.println("Sample: " + (i + 1));
				}
				if (frame != null) {
					String text = " Samples taken: " + (i+1);
					//remainingTime((currentTime - start)
					//		* ((period - i - 1) * sampRate
					//				+ sampRate - j - 1)
					//				/ (burnIn + i * sampRate + j + 1))
					if(autoPar.automateNumberOfSamplesToTake)
						text += " -- total number will be determined automatically";
					else
						text += " out of "+period;

					text += "   Sampling rate: " + sampRate;
					if(autoPar.automateNumberOfSamplesToTake){
						text +=  ",  Similarity(alignment n-1, alignment n): " + df.format(currScore) + " < " + df.format(AutomateParameters.PERCENT_CONST);
					}
					frame.statusText.setText(text );
				}
				if(autoPar.automateNumberOfSamplesToTake){
					alignment = new AlignmentData(getState().getLeafAlign());
					allAlignments.add(alignment);
					if (allAlignments.size() >1){
						FuzzyAlignment Fa = FuzzyAlignment.getFuzzyAlignmentAndProject(allAlignments.subList(0, allAlignments.size()-1), 0);
						FuzzyAlignment Fb = FuzzyAlignment.getFuzzyAlignmentAndProject(allAlignments, 0);
						//System.out.println(Fa);
						//System.out.println("xxxx");
						//System.out.println(Fb);
						currScore = FuzzyAlignment.AMA(Fa, Fb);
						System.out.println(currScore);
						distances.add(currScore);
						if (allAlignments.size() >5){
							shouldStop = AutomateParameters.shouldStopSampling(distances);
						}

					}
				}
				// Report the results of the sample.
				report(i, period);
			}
		} catch (StoppedException ex) {
			// stopped: report and save state
			errorCode = 1;
		}

		if(Utils.DEBUG) {
			System.out.println("Times spent in each MCMC step type (ms):");
			System.out.println(ali);
			System.out.println(top);
			System.out.println(edge);
			System.out.println(ind);
			System.out.println(sub);
		}

		// Triggers a /after first sample/ of the plugins.
		if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
			postprocMan.afterLastSample();
		}
		
		return errorCode;
	}

	private void doSwap(int swapNo) {
		int swapA, swapB;
		swapA = swapGenerator.nextInt(noOfProcesses);
		do {
			swapB = swapGenerator.nextInt(noOfProcesses);
		} while (swapA == swapB);

		System.out.printf("SwapNo: %d - SwapA: %d - SwapB: %d\n", swapNo,
				swapA, swapB);

		double swapAccept = swapGenerator.nextDouble();

		if (rank == swapA || rank == swapB) {
			double[] myStateInfo = new double[3];
			myStateInfo[0] = getTree().getLogLike();
			myStateInfo[1] = getTree().getLogPrior();
			myStateInfo[2] = heat;

			double[] partnerStateInfo = new double[3];

			mpi.Request send, recieve;

			if (rank == swapA) {
				send = MPI.COMM_WORLD.Isend(myStateInfo, 0, 3, MPI.DOUBLE,
						swapB, 0);
				recieve = MPI.COMM_WORLD.Irecv(partnerStateInfo, 0, 3,
						MPI.DOUBLE, swapB, 1);
			} else {
				send = MPI.COMM_WORLD.Isend(myStateInfo, 0, 3, MPI.DOUBLE,
						swapA, 1);
				recieve = MPI.COMM_WORLD.Irecv(partnerStateInfo, 0, 3,
						MPI.DOUBLE, swapA, 0);
			}

			mpi.Request.Waitall(new mpi.Request[] { send, recieve });

			System.out
			.printf("[Worker %d] Heat: [%f] - Sent: [%f,%f,%f] - Recv: [%f,%f,%f]\n",
					rank, heat, myStateInfo[0], myStateInfo[1],
					myStateInfo[2], partnerStateInfo[0],
					partnerStateInfo[1], partnerStateInfo[2]);

			double myLogLike = myStateInfo[0];
			double myLogPrior = myStateInfo[1];
			double myTemp = myStateInfo[2];
			double hisLogLike = partnerStateInfo[0];
			double hisLogPrior = partnerStateInfo[1];
			double hisTemp = partnerStateInfo[2];

			double acceptance = myTemp * (hisLogLike + hisLogPrior) + hisTemp
					* (myLogLike + myLogPrior);
			acceptance -= hisTemp * (hisLogLike + hisLogPrior) + myTemp
					* (myLogLike + myLogPrior);

			MPIUtils.println(rank,
					"Math.log(swapAccept): " + Math.log(swapAccept));
			MPIUtils.println(rank, "acceptance:           "
					+ acceptance);

			if (acceptance > Math.log(swapAccept)) {
				MPIUtils.println(rank,
						"Just swapped heat with my partner. New heat: "
								+ hisTemp);
				heat = hisTemp;
			}

			// MPI.COMM_WORLD.Send(myStateInfo, 0, 3, MPI.DOUBLE,
			// swapB, 0);
			// statalign.Utils.printLine(swapA, "Just sent " + swapB
			// + " my state.");
		}

	}

	private static String remainingTime(long x) {
		x /= 1000;
		return String.format("Estimated time left: %d:%02d:%02d", x / 3600,
				(x / 60) % 60, x % 60);
	}

	SimpleStats ali;
	SimpleStats top;
	SimpleStats edge;
	SimpleStats ind;
	SimpleStats sub;

	{
		if(Utils.DEBUG) {
			ali = new SimpleStats("Alignment");
			top = new SimpleStats("Topology");
			edge = new SimpleStats("Edge len");
			ind = new SimpleStats("Indel param");
			sub = new SimpleStats("Subst param");
		}
	}

	private void sample(int samplingMethod) throws StoppedException {
		if (samplingMethod == 0) {
			long timer;
			stoppable();
			switch (getTree().getSubstitutionModel().params != null
					&& getTree().getSubstitutionModel().params.length > 0 ? Utils
							.weightedChoose(FIVECHOOSE) : Utils
							.weightedChoose(FOURCHOOSE)) {
							case 0:
								if(Utils.DEBUG)
									timer = -System.currentTimeMillis();
                                edgeSampled++;
								edgeAccepted += strategy.sampleEdge() ? 1 : 0;
								if(Utils.DEBUG) {
									timer += System.currentTimeMillis();
									edge.addData(timer);
								}
								break;
							case 1:
								if(Utils.DEBUG)
									timer = -System.currentTimeMillis();
                                topologySampled++;
								topologyAccepted += strategy.sampleTopology() ? 1 : 0;
								if(Utils.DEBUG) {
									timer += System.currentTimeMillis();
									top.addData(timer);
								}
								break;
							case 2:
								if(Utils.DEBUG)
									timer = -System.currentTimeMillis();
                                indelSampled++;
								indelAccepted += strategy.sampleIndelParameter() ? 1 : 0;
								if(Utils.DEBUG) {
									timer += System.currentTimeMillis();
									ind.addData(timer);
								}
								break;
							case 3:
								if(Utils.DEBUG)
									timer = -System.currentTimeMillis();
                                alignmentSampled++;
								alignmentAccepted += strategy.sampleAlignment() ? 1 : 0;
								if(Utils.DEBUG) {
									timer += System.currentTimeMillis();
									ali.addData(timer);
								}
								break;
							case 4:
								if(Utils.DEBUG)
									timer = -System.currentTimeMillis();
                                substSampled++;
								substAccepted += strategy.sampleSubstParameter() ? 1 : 0;
								if(Utils.DEBUG) {
									timer += System.currentTimeMillis();
									sub.addData(timer);
								}
								break;
			}
		} else {
			stoppable();
			strategy.sampleEdge();
            strategy.sampleTopology();
            strategy.sampleIndelParameter();
            strategy.sampleSubstParameter();
            strategy.sampleAlignment();
		}
	}

	/**
	 * Returns a string representation describing the acceptance ratios of the current MCMC run.
	 * @return a string describing the acceptance ratios.
	 */
	public String getInfoString() {
		return String.format(Locale.US, "Acceptances: [Alignment: %f, Edge: %f, Topology: %f, Indel: %f, Substitution: %f]",
				(alignmentSampled == 0 ? 0 : (double) alignmentAccepted
						/ (double) alignmentSampled),
						(edgeSampled == 0 ? 0 : (double) edgeAccepted
								/ (double) edgeSampled),
								(topologySampled == 0 ? 0 : (double) topologyAccepted
										/ (double) topologySampled),
										(indelSampled == 0 ? 0 : (double) indelAccepted
												/ (double) indelSampled),
												(substSampled == 0 ? 0 : (double) substAccepted
														/ (double) substSampled));
	}

	/**
	 * Returns a {@link State} object that describes the current state of the
	 * MCMC. This can then be passed on to other classes such as postprocessing
	 * plugins.
	 */
	public State getState() {
		return getTree().getState();
	}

	private boolean isColdChain() {
		return heat == 1.0d;
	}

	private State MPIStateReceieve(int peer) {
        int tag = 0;

        int nn = 0, nl = 0;
        MPI.COMM_WORLD.Recv(nn, 0, 1, MPI.INT, peer, tag++);
        MPI.COMM_WORLD.Recv(nl, 0, 1, MPI.INT, peer, tag++);

        // Creates a new, uninitialized state and initializes the variables.
        State state = new State(nn, nl);

        // Set the names
        int[] nameLengths = new int[nn];
        for (int i = 0; i < nn; i++) {
            char[] name = new char[nameLengths[i]];
            MPI.COMM_WORLD.Recv(name, 0, nameLengths[i], MPI.CHAR, peer, tag++);
            state.name[i] = new String(name);
        }

        for (int i = 0; i < nn; i++)
            MPI.COMM_WORLD.Recv(state.children[i], 0, state.children[i].length, MPI.INT, peer, tag++);

		// parent
		MPI.COMM_WORLD.Recv(state.parent, 0, nn, MPI.INT, peer, tag++);
		// edgeLen
		MPI.COMM_WORLD.Recv(state.edgeLen, 0, nn, MPI.DOUBLE, peer, tag++);

		// sequences
		int[] seqLengths = new int[nn];
		MPI.COMM_WORLD.Recv(seqLengths, 0, nn, MPI.INT, peer, tag++);

		for (int i = 0; i < nn; i++) {
			char[] c = new char[seqLengths[i]];
			MPI.COMM_WORLD.Recv(c, 0, seqLengths[i], MPI.CHAR, peer, tag++);
			state.seq[i] = new String(c);
		}

		// align
		Object[] recvObj = new Object[1];
		MPI.COMM_WORLD.Recv(recvObj, 0, 1, MPI.OBJECT, peer, tag++);
		state.align = (int[][]) recvObj[0];

		// felsen
		MPI.COMM_WORLD.Recv(recvObj, 0, 1, MPI.OBJECT, peer, tag++);
		state.felsen = (double[][][]) recvObj[0];

		// indelParams
		final int noOfIndelParameter = 3;
		state.indelParams = new double[noOfIndelParameter];
		MPI.COMM_WORLD.Recv(state.indelParams, 0, noOfIndelParameter,
				MPI.DOUBLE, peer, tag++);

		// substParams
		int l = getTree().getSubstitutionModel().params.length;
		state.substParams = new double[l];
		MPI.COMM_WORLD.Recv(state.substParams, 0, l, MPI.DOUBLE, peer, tag++);

		// log-likelihood
		double[] d = new double[1];
		MPI.COMM_WORLD.Recv(d, 0, 1, MPI.DOUBLE, peer, tag++);
		state.logLike = d[0];

		// root
		int[] root = new int[1];
		MPI.COMM_WORLD.Recv(root, 0, 1, MPI.INT, peer, tag++);
		state.root = root[0];

		return state;
	}

	private void MPIStateSend(State state) {

		String[] seq = state.seq;
		int[][] align = state.align;
		double[][][] felsen = state.felsen;
		int nn = state.nn;
		int tag = 0;

        MPI.COMM_WORLD.Send(state.nn, 0, 1, MPI.INT, 0, tag++);
        MPI.COMM_WORLD.Send(state.nl, 0, 1, MPI.INT, 0, tag++);

        int[] nameLength = new int[nn];
        char[][] nameChars = new char[nn][];
        for (int i = 0; i < nn; i++) {
            nameLength[i] = state.name[i].length();
            nameChars[i] = state.name[i].toCharArray();
        }

        MPI.COMM_WORLD.Send(nameLength, 0, nn, MPI.INT, 0, tag++);
        for (int i = 0; i < nn; i++) {
            MPI.COMM_WORLD.Send(nameChars, 0, nameChars.length, MPI.CHAR, 0, tag++);
        }

        // TODO: Consider this!
        for (int i = 0; i < nn; i++)
            MPI.COMM_WORLD.Send(state.children[i], 0, state.children[i].length, MPI.INT, 0, tag++);

		// parent
		MPI.COMM_WORLD.Send(state.parent, 0, nn, MPI.INT, 0, tag++);
		// edgeLen
		MPI.COMM_WORLD.Send(state.edgeLen, 0, nn, MPI.DOUBLE, 0, tag++);

		// TODO: START OF OPTIMIZATION.

		// sequences
		int[] seqLength = new int[nn];
		char[][] seqChars = new char[nn][];
		for (int i = 0; i < nn; i++) {
			seqLength[i] = seq[i].length();
			seqChars[i] = seq[i].toCharArray();
		}
		MPI.COMM_WORLD.Send(seqLength, 0, nn, MPI.INT, 0, tag++);
		for (int i = 0; i < nn; i++) {
			MPI.COMM_WORLD.Send(seqChars[i], 0, seqLength[i], MPI.CHAR, 0, tag++);
		}

		// align
		Object[] alignObj = new Object[1];
		alignObj[0] = align;
		MPI.COMM_WORLD.Send(alignObj, 0, 1, MPI.OBJECT, 0, tag++);
		/*
		 * int[] alignLength = new int[align.length]; for (int i = 0; i <
		 * seq.length; i++) { alignLength[i] = align[i].length; }
		 * MPI.COMM_WORLD.Send(alignLength, 0, nn, MPI.INT, 0, tag++); for (int
		 * i = 0; i < align.length; i++) { MPI.COMM_WORLD.Send(align[i], 0,
		 * alignLength[i], MPI.INT, 0, tag++); }
		 */

		// felsen
		Object[] felsenObj = new Object[] { felsen };
		MPI.COMM_WORLD.Send(felsenObj, 0, 1, MPI.OBJECT, 0, tag++);

		// indelParams
		MPI.COMM_WORLD.Send(state.indelParams, 0, 3, MPI.DOUBLE, 0, tag++);

		// substParams
		MPI.COMM_WORLD.Send(state.substParams, 0, state.substParams.length,
				MPI.DOUBLE, 0, tag++);

		// loglikelihood
		MPI.COMM_WORLD.Send(new double[] { state.logLike }, 0, 1, MPI.DOUBLE,
				0, tag++);

		// root
		MPI.COMM_WORLD.Send(new int[] { state.root }, 0, 1, MPI.INT, 0, tag++);

		// TODO: END OF OPTIMIZATION.

	}

	private void report(int no, int total) {

		int coldChainLocation = -1;

		if (isParallel) {
			// Get rank of cold chain.
			int[] ranks = new int[] { (isColdChain() ? rank : 0) };
			int[] coldChainLoc = new int[1];
			MPI.COMM_WORLD.Reduce(ranks, 0, coldChainLoc, 0, 1, MPI.INT, MPI.SUM, 0);
			coldChainLocation = coldChainLoc[0];

			// TODO: Remove - for debugging purposes
			if (MPIUtils.isMaster(rank)) {
				MPIUtils.println(rank, "Cold chain is at: " + coldChainLocation);
			}

			if (isColdChain() && MPIUtils.isMaster(rank)) {
				// Sample normally.
				postprocMan.newSample(getState(), no, total);
			} else if (isColdChain() && !MPIUtils.isMaster(rank)) {
				// Send state.
				State state = getState();
				MPIStateSend(state);
			} else if (!isColdChain() && MPIUtils.isMaster(rank)) {
				// Receive state.
				State state = MPIStateReceieve(coldChainLocation);
				postprocMan.newSample(state, no, total);
			}

		} else {
			postprocMan.newSample(getState(), no, total);
		}

		// Log the accept ratios/params to the (.log) file. TODO: move to a plugin.
		try {
			if ((isParallel && MPIUtils.isMaster(rank)) || !isParallel) {
				postprocMan.logFile.write(getInfoString() + "\n");
				postprocMan.logFile.write("Report\tLogLikelihood\t"
						// + (tree.root.orphanLogLike + tree.root.indelLogLike)
                        + getTree().getLogLike()
						+ "\tR\t" + getTree().getR() + "\tLamda\t"
						+ getTree().getLambda() + "\tMu\t" + getTree().getMu()
								+ "\t" + getTree().getSubstitutionModel().print() + "\n");
				if (isParallel) {
					postprocMan.logFile.write("Cold chain location: " + coldChainLocation + "\n");
				}

			}
		} catch (IOException e) {
			if (postprocMan.mainManager.frame != null) {
				ErrorMessage.showPane(postprocMan.mainManager.frame, e, true);
			} else {
				e.printStackTrace(System.out);
			}
		}

		// alignmentSampled = 0;
		// alignmentAccepted = 0;
		// edgeSampled = 0;
		// edgeAccepted = 0;
		// topologySampled = 0;
		// topologyAccepted = 0;
		// indelSampled = 0;
		// indelAccepted = 0;
		// substSampled = 0;
		// substAccepted = 0;

	}

    public ITree getTree() {
        return strategy.getTree();
    }



	/**
	 * This function is only for testing and debugging purposes.
	 * 
	 * @param args
	 *            Actually, we do not use these parameters, as this function is
	 *            for testing and debugging purposes. All necessary input data
	 *            is directly written into the function.
	 * 
	 */
	// public static void main(String[] args) {
	// try {
	// Tree tree = new Tree(new String[] { "kkkkkkwwwwwwwwlidwwwwwkkk",
	// "kkkwwwwwwwlidwwwwwkkk", "kkkwwwwwwwlidwwwwwkkk",
	// "kkkwwwwwwwlidwwwwwkkk", "kkkwwwwwlidwwwwwkkkddkldkl",
	// "kkkwwwwwlidwwwwwkkkeqiqii", "kkkwwwwwlidwwwwwkkkddkidkil",
	// "kkkwwwwwlidwwwwwkkkeqiq", "kkkwwwwwlidwwwwwkkkddkldkll",
	// "kkkwwwwwlidwwwwwkkkddkldkil" }, new String[] { "A", "B",
	// "C", "D", "E", "F", "G", "H", "I", "J" }, new Dayhoff(), new Blosum62(),
	// "");
	// for (int i = 0; i < tree.vertex.length; i++) {
	// // tree.vertex[i].edgeLength = 0.1;
	// tree.vertex[i].updateTransitionMatrix();
	// }
	// tree.root.calcFelsRecursively();
	// System.out.println(tree.printedTree());
	// // Mcmc mcmc = new Mcmc(tree, new MCMCPars(0, 10000, 10, 1L), new
	// // PostprocessManager(null));
	// // mcmc.doMCMC();
	// } catch (StoppedException e) {
	// // stopped during tree construction
	// } catch (IOException e) {
	// }
	// }

}
