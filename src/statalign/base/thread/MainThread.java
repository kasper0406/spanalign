package statalign.base.thread;

import statalign.base.*;
import statalign.io.RawSequences;
import statalign.model.score.SubstitutionScore;
import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Kimura3;

import java.util.HashMap;
import java.util.Map;

/**
 * The main (suspendable) thread for background MCMC calculation.
 *
 * @author novak
 */
public class MainThread extends StoppableThread {
	
	/**
	 * Reference to the (singleton) MainManager object encapsulating all settings and data
	 * that an MCMC run depends on.
	 */
	public MainManager owner;
	
	/**
	 * Constructs a new MainThread that can be used to fire a background MCMC calculation.
	 * 
	 * @param owner Reference to the MainManager object.
	 */
	public MainThread(MainManager owner) {
		this.owner = owner;
	}
	/**
	 * Start background MCMC calculation.
	 */
	@Override
	public synchronized void run() {
		try {

			owner.postProcMan.initRun(owner.inputData);
			
			RawSequences seqs = owner.inputData.seqs;
			
			if(owner.frame != null) {
				owner.frame.statusText.setText(" Generating initial tree and alignment...");
			}

			System.out.println("\nPreparing initial tree and alignment...\n");

			// remove gaps and whitespace
			String[] nongapped = new String[seqs.size()];
			StringBuilder builder = new StringBuilder();
			int i, j;
			char ch;
			for(i = 0; i < nongapped.length; i++) {
				builder.setLength(0);
				String seq = seqs.getSequence(i);
				for(j = 0; j < seq.length(); j++) {
					ch = seq.charAt(j);
					if(Character.isWhitespace(ch) || ch == '-')
						continue;
					builder.append(ch);
				}
				nongapped[i] = builder.toString();
			}

            Mcmc mcmc;
            switch (owner.inputData.pars.treeType) {
                case STEINER:
                    Tree tree = new Tree(nongapped, seqs.getSeqnames().toArray(new String[seqs.size()]),
                            owner.inputData.model,
                            owner.inputData.model.attachedScoringScheme);
                    mcmc = new Mcmc(new SteinerTreeMCMCStrategy(tree), owner.inputData.pars, owner.postProcMan);
                    break;

                case SPANNOID:
                    int componentSize = owner.inputData.pars.componentSize;
                    Spannoid.BonphyStrategy bonphyStrategy = owner.inputData.pars.bonphyStrategy;
                    Spannoid spannoid = new Spannoid(componentSize, bonphyStrategy, nongapped,
                            seqs.getSeqnames().toArray(new String[seqs.size()]),
                            owner.inputData.model, owner.inputData.model.attachedScoringScheme);
                    mcmc = new Mcmc(new SpannoidMCMCStrategy(spannoid, owner.inputData.pars.componentSize), owner.inputData.pars, owner.postProcMan);
                    break;

                default:
                    throw new RuntimeException("Invalid tree topology!");
            };

            /*
            String newick = "((D:0.1,(F:0.1,(E:0.2,C:0.05):0.2)B:0.1):0.2)A;";
            String[] my_seqs = new String[] { "AAGT", "CGATTC", "CCGAAG", "AGACA", "TTGACC", "GTAC" };
            Map<String, Integer> nameMap = new HashMap<String, Integer>();
            for (int k = 0; k < my_seqs.length; k++)
                nameMap.put("" + (char)((int)'A' + k), k);

            ITree tree = new Spannoid(newick, my_seqs, nameMap,
                    owner.inputData.model, owner.inputData.model.attachedScoringScheme);
                    */

			int errorCode = mcmc.doMCMC();

            /*
            // TODO: Remove this
            //       Hack to print structure before MCMC.
            owner.postProcMan.beforeFirstSample();
            // owner.postProcMan.newStep(new McmcStep());
            owner.postProcMan.newSample(tree.getState(), 0, 0);
            owner.postProcMan.afterLastSample();
            */

			owner.postProcMan.finalizeRun();
			owner.finished(errorCode, null);
			// System.out.println(errorCode == 0 ? "Ready." : "Stopped.");
		} catch (StoppedException e) {
			// stopped during initial alignment
			owner.postProcMan.finalizeRun();
			owner.finished(2, null);
			System.out.println("Stopped.");
		} catch(Exception e) {
			owner.postProcMan.finalizeRun();
			e.printStackTrace();
			owner.finished(-1, e);
		}
		
	}
}
