package statalign.postprocess.plugins;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import statalign.base.CircularArray;
import statalign.base.InputData;
import statalign.base.State;
import statalign.base.Utils;
import statalign.postprocess.Postprocess;
import statalign.postprocess.gui.AlignmentGUI;

public class MpdAlignment extends statalign.postprocess.Postprocess {

	public String title;	
	public int frequency = 5;
	JPanel pan;
	AlignmentGUI gui;
    AlignmentGUI gui_ref;
	//private boolean sampling = true;

	// CurrentAlignment curAlig;

	ColumnNetwork network;
	Column firstVector, lastVector;
	int sizeOfAlignments;

	int[] firstDescriptor; 
	String t[][];
	String[] sequences;
	String[] viterbialignment;

	double[] decoding;
	String[] alignment;

	InputData input;

    AlignmentScorer scorer;

	public MpdAlignment(){
		screenable = true;
		outputable = true;
		postprocessable = true;
		postprocessWrite = true;
		rnaAssociated = false;
	}

	@Override
	public JPanel getJPanel() {
		pan = new JPanel(new GridLayout(0,1));
		return pan;
	}

	@Override
	public Icon getIcon() {
//		return new ImageIcon("icons/MPD.gif");
		return new ImageIcon(ClassLoader.getSystemResource("icons/MPD.gif"));
	}

	@Override
	public String getTabName() {
		return "MPD";
	}

	@Override
	public String getTip() {
		return "Maximum Posterior Decoding (consensus) alignment";
	}

	@Override
	public double getTabOrder() {
		return 5.0d;
	}

	@Override
	public String getFileExtension() {
		return "mpd";
	}

	@Override
	public String[] getDependences() {
		return new String[] { "statalign.postprocess.plugins.CurrentAlignment" };
	}


	@Override
	public void refToDependences(Postprocess[] plugins) {
		// curAlig = (CurrentAlignment) plugins[0];
	}

	static Comparator<String[]> compStringArr = new Comparator<String[]>() {
		@Override
		public int compare(String[] a1, String[] a2) {
			return a1[0].compareTo(a2[0]);
		}};

	@Override
	public void beforeFirstSample(InputData input) {
		if(show) {
			pan.removeAll();
			title = input.title;
			JScrollPane scroll = new JScrollPane();
			scroll.setViewportView(gui = new AlignmentGUI(title,input.model));//, mcmc.tree.printedAlignment()));
			pan.add(scroll);
            scroll = new JScrollPane();
            scroll.setViewportView(gui_ref = new AlignmentGUI("Reference alignment",input.model));
            pan.add(scroll);
//			System.out.println("Mpd Alignment parent: " + pan.getParent());
			pan.getParent().validate();
		}

		this.input = input;

		sizeOfAlignments = input.seqs.size();
		alignment = new String[sizeOfAlignments];
		if(show)
			gui.alignment = alignment;
		t = new String[sizeOfAlignments][];
		sequences = null;
		viterbialignment = new String[sizeOfAlignments];

        String[] sorted = new String[input.seqs.size()];
        for (int i = 0; i < input.seqs.size(); i++) {
            sorted[i] = input.seqs.getSeqNamePadded(i) + "\t" + input.seqs.getSequence(i);
        }
        Arrays.sort(sorted);

        String[][] temp = new String[sizeOfAlignments][];
        for (int i = 0; i < sizeOfAlignments; i++){
            temp[i] = sorted[i].split("\t");
        }
        Arrays.sort(temp, compStringArr);
        String[] ref = new String[sizeOfAlignments];
        for (int i = 0; i < sizeOfAlignments; i++){
            ref[i] = temp[i][1];
        }
        scorer = new AlignmentScorer(ref);

        if (show)
            gui_ref.alignment = sorted;

		network = new ColumnNetwork();

		firstDescriptor = new int[sizeOfAlignments];
		Arrays.fill(firstDescriptor, -1);
		firstVector = network.add(firstDescriptor);
		lastVector = null;
	}

	@Override
	public void newSample(State state, int no, int total) {
		//System.out.println(curAlig);
		//System.out.println(curAlig.leafAlignment);
        String[] leafAlignment = new String[state.nl];
        String[] rows = state.getLeafAlign();
        String[] names = state.getPaddedNames();
        for (int i = 0; i < rows.length; i++)
            leafAlignment[i] = names[i] + '\t' + rows[i];

		for(int i = 0; i < t.length; i++){
            if (leafAlignment[i] == null) {
                System.out.println();
            }
			t[i] = leafAlignment[i].split("\t");
		}
		Arrays.sort(t, compStringArr);

		int[] previousDescriptor = firstDescriptor;

		int i, j, len = t[0][1].length();
		for(j = 0; j < len; j++){
			int[] nextDescriptor =  new int[sizeOfAlignments];
			boolean allGap = true;
			for(int k = 0; k < sizeOfAlignments; k++){
				if(t[k][1].charAt(j) == '-')
					nextDescriptor[k] = ColumnKey.colNext(previousDescriptor[k]);
				else {
					nextDescriptor[k] = ColumnKey.colNext(previousDescriptor[k])+1;
					allGap = false;
				}
			}
			if(!allGap)
				network.add(nextDescriptor);//[j]);

			previousDescriptor = nextDescriptor;
		}//j (length of alignments)

		if(no == 0) {		// add last vector once only
			int[] lastDescriptor =  new int[sizeOfAlignments];
			for(j = 0; j < sizeOfAlignments; j++){
				lastDescriptor[j] = ColumnKey.colNext(previousDescriptor[j])+1;
			}
			lastVector = network.add(lastDescriptor);
		}
		if(no == 0 || (total-1-no) % frequency == 0) {
			network.updateViterbi(no+1);
			//System.out.println("sequences first: "+sequences);

			if(sequences == null) {
				sequences = new String[sizeOfAlignments];
				for(i = 0; i < sizeOfAlignments; i++){
					sequences[i] = "";
					for(j = 0; j < len; j++){
						if(t[i][1].charAt(j) != '-'){
							sequences[i] += t[i][1].charAt(j);
						}
					}
				}
			}
			for(i = 0; i < sizeOfAlignments; i++)
				viterbialignment[i] = "";
			Column actualVector = lastVector.viterbi;
			ArrayList<Integer> posteriorList = new ArrayList<Integer>();
			while(!actualVector.equals(firstVector)){
				int[] desc = actualVector.key.desc;
				posteriorList.add(new Integer(actualVector.count));
				for(i = 0; i < desc.length; i++){
					if((desc[i] & 1) == 0){
						viterbialignment[i] = "-"+viterbialignment[i];
					}
					else{
						viterbialignment[i] = sequences[i].charAt(desc[i] >> 1) + viterbialignment[i];
					}
				}
				actualVector = actualVector.viterbi;
			}
			decoding = new double[posteriorList.size()];
			for(i = 0; i < decoding.length; i++){
				decoding[i] = (double)(posteriorList.get(posteriorList.size() - i - 1)).intValue()/(no+1);
			}

			for(i = 0; i < viterbialignment.length; i++){
				alignment[i] = t[i][0]+"\t"+viterbialignment[i];
			}

			// sort alignment lexicographically
			// TODO sort oder is parameter (alternatives: original, tree, lexico)
			Arrays.sort(alignment);

			if(show) {
				gui.decoding = decoding;
				gui.alignment = alignment;
				gui.repaint();

                gui_ref.title = String.format("Score against reference alignment: %.3f", scorer.calculateScore(viterbialignment));
                gui_ref.repaint();
			}
		}
		
		if(sampling){
			try {
				String[] aln = Utils.alignmentTransformation(alignment, alignmentType, input);
				for(i = 0; i < aln.length; i++){
					file.write("Sample "+no+"\tMPD alignment:\t"+aln[i]+"\n");
				}
				if(decoding != null){
					for(i = 0; i < decoding.length; i++){
						file.write("Sample "+no+"\tMPD alignment probabilities:\t"+decoding[i]+"\n");
					}
				}
				else{
					file.write("Sample "+no+"\tMPD alignment:\tNo posterior values so far\n");
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}	



		/*if(no == 0)
		{
			String [] [] inputAlignment = new String[input.seqs.sequences.size()][2];
			for(int k = 0 ; k < inputAlignment.length ; k++)
			{
				inputAlignment[k][0] = input.seqs.seqNames.get(k);
				inputAlignment[k][1] = input.seqs.sequences.get(k);
			}
			Arrays.sort(inputAlignment, compStringArr);
			
			appendAlignment("reference", inputAlignment, new File(input.title+".samples"), false);
			appendAlignment(no+"", t, new File(input.title+".samples"), true);
		}
		else
		{
			appendAlignment(no+"", t, new File(input.title+".samples"), true);
		}*/

	}

    public double getScore() {
        return scorer.calculateScore(viterbialignment);
    }

    @Override
	public void afterLastSample() {
		if (postprocessWrite) {
			try {
				String[] aln = Utils.alignmentTransformation(alignment,
						alignmentType, input);
				for (int i = 0; i < aln.length; i++) {
					outputFile.write(aln[i] + "\n");
				}
				outputFile.write("\n#scores\n\n");
				if (decoding != null) {
					for (int i = 0; i < decoding.length; i++) {
						outputFile.write(decoding[i] + "\n");
					}
				} else {
					outputFile.write("No posterior values so far\n");
				}
			} catch (IOException e) {
			}
			/*
			appendAlignment("mpd", alignment, new File(input.title+".samples"), true);
			PPFold.saveToFile(Utils.alignmentTransformation(alignment,
					"Fasta", input), new File(input.title+".dat.res.mpd"));
			try
			{

				BufferedWriter buffer = new BufferedWriter(new FileWriter(new File(input.title+".samples"), true));
				buffer.write("%posteriors\n");
				for(int i = 0 ; i < decoding.length ; i++)
				{
					buffer.write(decoding[i]+"\n");
				}
				buffer.close();
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}*/
		}
	}
	
	@Override
	public List<String> getCreatedFileDescriptions() {
		return Arrays.asList("The MPD (consensus) alignment");
	}

	/* (non-Javadoc)
	 * @see statalign.postprocess.Postprocess#setSampling(boolean)
	 */
	@Override
	public void setSampling(boolean enabled) {
		sampling = enabled;

	}

}

class ColumnNetwork {
	HashMap<ColumnKey,Column> contMap = new HashMap<ColumnKey,Column>();
	HashMap<ColumnKey,ArrayList<Column>> preMap = new HashMap<ColumnKey,ArrayList<Column>>();
	HashMap<ColumnKey,ArrayList<Column>> postMap = new HashMap<ColumnKey,ArrayList<Column>>();

	int numberOfEdges = 0;
	int numberOfNodes = 0;

	Column first;

	/**
	 * Adds a new alignment column into the network. If already in the network, MyVector.count is incremented.
	 * @param descriptor Alignment column represented by an array of signed integers
	 * @return MyVector of alignment column or null if it was in network before
	 */
	Column add(int[] descriptor) {
		Column val;

		ColumnKey key = new ColumnKey(descriptor);
		if((val=contMap.get(key)) != null) {
			val.logcnt = Math.log(++val.count);
			return null;
		}
		val = new Column(key);
		contMap.put(key, val);
		if(numberOfNodes == 0)
			first = val;
		numberOfNodes++;

		ArrayList<Column> arr;

		key = new ColumnKey(ColumnKey.pre(descriptor));
		if((arr = preMap.get(key)) == null)
			preMap.put(key, arr = new ArrayList<Column>());
		arr.add(val);
		if((arr = postMap.get(key)) != null) {
			for(Column postCol : arr) {
				postCol.inNum++;
				val.outgoing.add(postCol);
				numberOfEdges++;
			}
		}

		key = new ColumnKey(ColumnKey.post(descriptor));
		if((arr = postMap.get(key)) == null)
			postMap.put(key, arr = new ArrayList<Column>());
		arr.add(val);
		if((arr = preMap.get(key)) != null) {
			for(Column preCol : arr) {
				preCol.outgoing.add(val);
				val.inNum++;
				numberOfEdges++;
			}
		}

		return val;
	}


	void updateViterbi(int n) {
		double logN = Math.log(n);
		CircularArray<Column> calculable = new CircularArray<Column>();
		calculable.push(first);
		first.score = logN;
		Column act;
		while((act = calculable.shift()) != null) {
			double myScore = act.score+act.logcnt-logN;
			for(Column outGoing : act.outgoing) {
				if(outGoing.score < myScore) {
					outGoing.score = myScore;
					outGoing.viterbi = act;
				}
				if(++outGoing.inReady == outGoing.inNum)
					calculable.push(outGoing);
			}
			act.inReady = 0;
			act.score = -1e300;
		}
	}

}

class Column {
	ArrayList<Column> outgoing = new ArrayList<Column>();
	ColumnKey key;
	int inNum;

	int inReady;
	int count = 1;
	double logcnt = 0;
	double score = -1e300;
	Column viterbi;

	Column(ColumnKey _key) {
		key = _key;
	}

}

class ColumnKey {
	public int[] desc;

	ColumnKey(int[] arr) {
		desc = arr;
	}

	@Override
	public boolean equals(Object o)	{
		return (o instanceof ColumnKey) && Arrays.equals(desc, ((ColumnKey)o).desc);
	}

	@Override
	public int hashCode()	{
		return Arrays.hashCode(desc);
	}

	static int[] pre(int[] desc) {
		int[] ret = new int[desc.length];
		for(int i = 0; i < desc.length; i++)
			ret[i] = (desc[i]+1) >> 1;
		return ret;
	}

	static int[] post(int[] desc) {
		int[] ret = new int[desc.length];
		for(int i = 0; i < desc.length; i++)
			ret[i] = desc[i] >> 1;
		return ret;
	}

	static int colNext(int n) {
		return n + (n & 1);
	}

}

class AlignmentScorer {
    private String[] ref;
    int[] lengths;

    AlignmentScorer(String[] ref) {
        this.ref = ref;

        lengths = new int[ref.length];
        for (int i = 0; i < ref.length; i++) {
            lengths[i] = ref[i].replace("-","").length();
        }
    }

    public double calculateScore(String[] pred) {
        int score = 0;
        int total = 0;

        for (int m = 0; m < ref.length-1; m++) {
            for (int n = m+1; n < ref.length; n++) {
                score += calculateScore(ref[m], ref[n], pred[m], pred[n]);
                total += lengths[m];

                score += calculateScore(ref[n], ref[m], pred[n], pred[m]);
                total += lengths[n];
            }
        }

        return (double)score / total;
    }

    private static int calculateScore(String refA, String refB, String predA, String predB) {
        List<Integer> ref = transcode(refA, refB);
        List<Integer> pred = transcode(predA, predB);

        int score = 0;
        for (int i = 0; i < ref.size(); i++) {
            if (ref.get(i).equals(pred.get(i))) {
                score++;
            }
        }

        return score;
    }

    private static List<Integer> transcode(String A, String B) {
        List<Integer> result = new ArrayList<Integer>();
        int counter = 0;
        for (int i = 0; i < A.length(); i++) {
            boolean charAtA = A.charAt(i) != '-';
            boolean charAtB = B.charAt(i) != '-';

            if (charAtA) {
                result.add(charAtB ? counter : -1);
            }
            if (charAtB) {
                counter++;
            }
        }
        return result;
    }
}