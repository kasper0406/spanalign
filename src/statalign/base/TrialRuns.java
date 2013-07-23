package statalign.base;

/**
 * Created with IntelliJ IDEA.
 * User: aldo
 * Date: 23/07/13
 * Time: 10:35
 * To change this template use File | Settings | File Templates.
 */
public class TrialRuns {
    public static void main(String[] args) throws Exception {
        for(int i=0; i< 1000; i++){System.out.println(Utils.generator.nextInt(3));}
    }

}
